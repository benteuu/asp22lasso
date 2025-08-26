set.seed(969898)

n <- 1000
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
y <- rnorm(n, 0 + 1 * x1 + 1 * x3, exp(-3 + 1 * x2 + 1 * x3))
m <- lasso_lmls(y ~ x1 + x3,
                ~ x2 + x3,
                data = c(y, x1, x2, x3))

# score -----------------------------------------------------------------------

beta <- gamma <- c(0, 0, 0)
m <- set_coef(m, "location", beta)
m <- set_coef(m, "scale", gamma)

f <- function(x, predictor) {
  loglik(set_coef(m, predictor, x))
}

test_that("score of beta is correct", {
  num_score <- numDeriv::grad(f, beta, predictor = "location")
  expect_exactly(score(m, "location"), num_score)
})

test_that("score of gamma is correct", {
  num_score <- numDeriv::grad(f, gamma, predictor = "scale")
  expect_exactly(score(m, "scale"), num_score)
})
