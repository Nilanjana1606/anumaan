#!/usr/bin/env Rscript
source_helpers <- function() {
 tail_probs <- function(obs,reps){s<-length(reps);lo<-(1+sum(reps<=obs))/(s+1);hi<-(1+sum(reps>=obs))/(s+1);c(lower_tail=lo,upper_tail=hi,p_two_sided=min(1,2*min(lo,hi)))}
 auc <- function(y,p){if(length(unique(y))<2L)return(NA_real_);r<-rank(p,ties.method="average");(sum(r[y==1])-sum(y==1)*(sum(y==1)+1)/2)/(sum(y==1)*sum(y==0))}
 list(tail_probs=tail_probs,auc=auc)
}
 h <- source_helpers(); stopifnot(abs(h$tail_probs(3,1:5)["lower_tail"]-4/6)<1e-12, h$tail_probs(0,1:5)["lower_tail"]==1/6, h$tail_probs(6,1:5)["upper_tail"]==1/6, h$tail_probs(3,rep(3,5))["p_two_sided"]==1, is.na(h$auc(c(1,1),c(.2,.8))), abs(h$auc(c(0,1),c(.1,.9))-1)<1e-12)
X<-diag(3); mask<-matrix(c(TRUE,FALSE,TRUE,TRUE,TRUE,FALSE),3,2); stopifnot(sum(mask[,1])==2, sum(mask[,2])==2, all(c("lower_tail","upper_tail","p_two_sided") %in% names(h$tail_probs(.5,c(.2,.5,.8)))))
cat("Q200 failure-decomposition deterministic tests: PASS\n")
