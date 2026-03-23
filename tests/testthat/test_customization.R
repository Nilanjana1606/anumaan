# Test Plot Customization Parameters

test_that("plot_bar accepts order parameter and sets title", {
  test_data <- data.frame(
    organism = rep(c("E. coli", "K. pneumoniae", "S. aureus"), c(30, 20, 10))
  )

  p <- plot_bar(
    data = test_data,
    x = "organism",
    title = "Test bar plot",
    order = "desc"
  )
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Test bar plot")

  p_asc <- plot_bar(
    data = test_data,
    x = "organism",
    order = "asc"
  )
  expect_s3_class(p_asc, "ggplot")
})

test_that("plot_bar with flip_coords works", {
  test_data <- data.frame(
    organism = rep(c("A", "B", "C"), c(50, 30, 10))
  )

  p <- plot_bar(
    data = test_data,
    x = "organism",
    flip_coords = TRUE,
    title = "Flipped",
    xlab = "Organism",
    ylab = "Count"
  )
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Flipped")
})

test_that("plot_proportion creates a plot", {
  test_data <- data.frame(
    organism = rep(c("E. coli", "K. pneumoniae"), each = 20),
    result = rep(c("R", "S"), 20)
  )

  p <- plot_proportion(
    data = test_data,
    x = "organism",
    fill = "result"
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_stacked_bar creates a plot", {
  test_data <- data.frame(
    drug_class = rep(c("Carbapenems", "Fluoroquinolones"), each = 20),
    result = rep(c("R", "S"), 20)
  )

  p <- plot_stacked_bar(
    data = test_data,
    x = "drug_class",
    fill = "result"
  )
  expect_s3_class(p, "ggplot")
})
