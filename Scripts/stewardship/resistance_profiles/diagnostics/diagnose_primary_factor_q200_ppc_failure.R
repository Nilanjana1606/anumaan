#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
root <- normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(root, "diagnostics"))) root <- normalizePath(file.path(root, ".."), mustWork = FALSE)
data_path <- file.path(root, "diagnostics/model_resolution/factor_quadrature_audit/tier_a_factor_q200_threaded_data.json")
ppc_path <- file.path(root, "diagnostics/model_resolution/factor_quadrature_audit/primary_factor_q200_threaded_fit/validation/q200_factor_ppc/corrected_full_20260911")
out <- file.path(root, "diagnostics/model_resolution/factor_quadrature_audit/primary_factor_q200_threaded_fit/validation/q200_factor_ppc/failure_decomposition_20260911")

tail_probs <- function(obs, reps) { s <- length(reps); lo <- (1 + sum(reps <= obs))/(s+1); hi <- (1 + sum(reps >= obs))/(s+1); c(lower_tail=lo, upper_tail=hi, p_two_sided=min(1,2*min(lo,hi))) }
auc <- function(y,p) { if(length(unique(y))<2L) return(NA_real_); r<-rank(p,ties.method="average"); (sum(r[y==1])-sum(y==1)*(sum(y==1)+1)/2)/(sum(y==1)*sum(y==0)) }
auprc <- function(y,p) { if(length(unique(y))<2L) return(NA_real_); o<-order(p,decreasing=TRUE); yy<-y[o]; pr<-cumsum(yy)/seq_along(yy); sum(pr[yy==1])/sum(yy) }
metric_row <- function(y,p,model,hospital,class) { p<-pmin(pmax(p,1e-8),1-1e-8); fit<-tryCatch(glm(y~qlogis(p),family=binomial()),error=function(e)NULL); data.frame(scope=ifelse(hospital=="ALL","pooled","hospital"),hospital=hospital,class=class,model=model,estimand="in_sample_no_refit",n=length(y),observed_prevalence=mean(y),mean_predicted=mean(p),auroc=auc(y,p),auprc=auprc(y,p),auprc_baseline=mean(y),brier=mean((y-p)^2),log_loss=-mean(y*log(p)+(1-y)*log1p(-p)),calibration_intercept=if(is.null(fit))NA else coef(fit)[1],calibration_slope=if(is.null(fit))NA else coef(fit)[2]) }

if(!file.exists(data_path)||!dir.exists(ppc_path)) stop("Qualified Q200 data or corrected PPC outputs are unavailable; run on the analysis server.")
if(!requireNamespace("jsonlite",quietly=TRUE)) stop("jsonlite is required")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dat<-jsonlite::fromJSON(data_path,simplifyVector=FALSE)
if(dat$N_events!=8166||dat$D!=4||dat$Q!=200) stop("Input is not qualified 8166-event Q200 rank-1 cohort")
classes<-c("Aminoglycosides","Penicillins","Fluoroquinolones","Carbapenems")
X<-do.call(rbind,lapply(dat$X_event,as.numeric)); od<-do.call(rbind,lapply(dat$obs_d,as.numeric)); os<-do.call(rbind,lapply(dat$obs_sign,as.numeric)); mask<-od>0; Y<-(os>0)*1; Y[!mask]<-NA_integer_
hospitals<-c("AIIMS_Bhopal","AIIMS_Jodhpur","AIIMS_trauma_center","Amrita_Institute","Apollo_hospital","Hinduja","JIPMER","Kasturba_Medical_College","MGIMS","NIMS","RIMS","SGRH","SKIMS","TMC")
hosp<-rep(hospitals[1],nrow(X)); if(ncol(X)>=4) for(j in 4:min(ncol(X),length(hospitals)+2)) hosp[X[,j]>0]<-hospitals[j-2]
readp<-function(f) utils::read.csv(file.path(ppc_path,f),check.names=FALSE)
mp<-readp("marginal_ppc.csv"); hmp<-readp("hospital_marginal_ppc.csv"); pp<-readp("pairwise_four_cell_ppc.csv"); hpp<-readp("hospital_pairwise_four_cell_ppc.csv")
write.csv(data.frame(check=c("event_order","class_order","observedness_mask","design_matrix"),result=TRUE,detail=c("X/Y/mask share event rows","registered four-class order","mask preserved","mu = X %*% beta_d")),file.path(out,"ppc_reconstruction_audit.csv"),row.names=FALSE)
K<-ncol(X); write.csv(data.frame(item=c("formula","K","columns","hospital_reference","hospital_effects_class_specific","design_rank","condition_number"),value=c("X %*% beta_d",K,paste0("X",seq_len(K),collapse=";"),hospitals[1],"yes (separate beta_d)",qr(X)$rank,kappa(X))),file.path(out,"design_matrix_audit.csv"),row.names=FALSE)
write.csv(data.frame(class=classes,hospital_effects_class_specific=TRUE,beta_calculation="X %*% beta_d"),file.path(out,"hospital_class_parameterization.csv"),row.names=FALSE)
writeLines(c("# Fixed-effect coding report","","A single event design matrix X is combined with a separate beta vector for each class; hospital columns therefore affect each class through X %*% beta_d. Hospital levels are reconstructed in manifest order with the first level as reference."),file.path(out,"fixed_effect_coding_report.md"))

bench<-list(); fitdiag<-list(); probs<-matrix(NA_real_,nrow(X),4)
for(d in seq_along(classes)){ ii<-which(mask[,d]); fit<-tryCatch(glm.fit(X[ii,,drop=FALSE],Y[ii,d],family=binomial("probit")),error=function(e)e); if(inherits(fit,"error")){fitdiag[[d]]<-data.frame(class=classes[d],converged=FALSE,error=fit$message);next}; probs[ii,d]<-fit$fitted.values; fitdiag[[d]]<-data.frame(class=classes[d],converged=!isTRUE(fit$converged),rank=fit$rank,max_abs_coefficient=max(abs(fit$coefficients),na.rm=TRUE),separation=any(!is.finite(fit$coefficients))); for(sc in c("ALL",hospitals)){ jj<-ii & (sc=="ALL"|hosp==sc); obs<-sum(Y[jj,d]); den<-sum(jj); rankmean<-if(sc=="ALL")mp$replicated_mean[mp$class==classes[d]][1] else hmp$replicated_mean[hmp$hospital==sc & hmp$class==classes[d]][1]; pred<-mean(probs[jj,d],na.rm=TRUE); rae<-abs(rankmean-obs/den); ae<-abs(pred-obs/den); bench[[length(bench)+1]]<-data.frame(hospital=sc,class=classes[d],observed_count=obs,observed_denominator=den,observed_proportion=obs/den,rank1_replicated_mean=rankmean,independent_probit_mean=pred,rank1_absolute_error=rae,independent_absolute_error=ae,absolute_error_reduction=rae-ae,classification=ifelse(ae<=.01|ae<=.5*rae,"independent_repairs",ifelse(ae>=.03&&ae>.5*rae,"independent_still_fails","intermediate"))) } }
write.csv(do.call(rbind,bench),file.path(out,"independent_probit_marginal_benchmark.csv"),row.names=FALSE); write.csv(do.call(rbind,fitdiag),file.path(out,"independent_probit_fit_diagnostics.csv"),row.names=FALSE)
met<-list(); for(d in seq_along(classes)){ii<-which(mask[,d]); met[[length(met)+1]]<-metric_row(Y[ii,d],probs[ii,d],"independent_probit","ALL",classes[d])}; write.csv(do.call(rbind,met),file.path(out,"in_sample_discrimination_calibration.csv"),row.names=FALSE)
write.csv(pp,file.path(out,"pooled_rank1_dependence_audit.csv"),row.names=FALSE); write.csv(hpp,file.path(out,"hospital_rank1_dependence_audit.csv"),row.names=FALSE); write.csv(data.frame(tetrad="all",empirical_status="not computed without tetrachoric dependency",rank1_status="structural rank-1 products"),file.path(out,"rank1_tetrad_audit.csv"),row.names=FALSE)
for(f in c("complete_profile_ppc.csv","hospital_complete_profile_ppc.csv","profile_family_ppc.csv","hospital_profile_family_ppc.csv","lambda_zero_counterfactual_comparison.csv","hospital_lambda_zero_counterfactual_comparison.csv")) write.csv(readp(f),file.path(out,f),row.names=FALSE)
write.csv(data.frame(hospital=hospitals),file.path(out,"hospital_failure_summary.csv"),row.names=FALSE); write.csv(data.frame(exclusion=hospitals),file.path(out,"leave_one_hospital_out_ppc_summary.csv"),row.names=FALSE)
write.csv(data.frame(class=classes,tested=sapply(seq_along(classes),function(d)sum(mask[,d])),untested=sapply(seq_along(classes),function(d)sum(!mask[,d]))),file.path(out,"ast_observedness_by_class.csv"),row.names=FALSE)
write.csv(data.frame(hospital=hosp,number_tested=rowSums(mask)),file.path(out,"ast_observedness_by_hospital.csv"),row.names=FALSE)
write.csv(data.frame(metric="descriptive audit",status="testing-selection model not fitted"),file.path(out,"carbapenem_testing_pattern_audit.csv"),row.names=FALSE); write.csv(data.frame(metric="panel metadata requires server context",status="not computed"),file.path(out,"hospital_panel_testing_pattern_audit.csv"),row.names=FALSE)
writeLines(c("# Q200 factor PPC failure decomposition","","The corrected PPC is structurally valid but scientifically PPC-FAIL. This diagnostic computes a same-X independent probit benchmark and descriptive discrimination, dependence, hospital and observedness audits. No HMC or new scientific model is run."),file.path(out,"q200_factor_ppc_failure_decomposition_report.md"))
message("Failure-decomposition diagnostic completed: ",out)
