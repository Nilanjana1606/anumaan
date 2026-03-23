# Test Plot Ordering Functionality

test_that("plot_bar orders bars by frequency descending by default", {
  test_data <- data.frame(
    organism = rep(c("A", "B", "C"), c(50, 30, 10))
  )

  p <- plot_bar(data = test_data, x = "organism", flip_coords = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_bar with flip_coords and title works", {
  test_data <- data.frame(
    organism = rep(c("E. coli", "K. pneumoniae", "S. aureus"), c(30, 20, 10))
  )

  p <- plot_bar(
    data = test_data,
    x = "organism",
    flip_coords = TRUE,
    title = "Top Organisms (Ordered by Frequency)",
    xlab = "Organism",
    ylab = "Count"
  )
  expect_s3_class(p, "ggplot")
})
