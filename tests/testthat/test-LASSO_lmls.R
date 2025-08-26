set.seed(969898)

n <- 1000
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
y <- rnorm(n, 0 + 1 * x1 + 1 * x3, exp(-3 + 1 * x2 + 1 * x3))
m <- lasso_lmls(y ~ x1 + x3,
                ~ x2 + x3,
                data = c(y, x1, x2, x3))

# Integration Tests:

# Test if setup works correctly -----------------------------------------------
m <- setup(y ~ x1 + x3,
           ~ x2 + x3,
           const = 1e-06)

test_that("setup() creates list with correct names", {
  n <- c("y", "x", "z", "lasso", "lambda", "const", "nobs", "df", "df.residual",
         "light", "call", "term", "coefficients", "fitted.values", "residuals",
         "vcov", "chol_info_gamma", "chol_info_beta", "iterations", "loss")

  expect_true(all(n %in% names(m)))
})

test_that("setup() creates list with correct values", {
  expect_equal(m$y, y)

  expect_exactly(m$x, cbind(1, x1, x3))
  expect_exactly(m$z, cbind(1, x2, x3))

  expect_equal(m$nobs, n)

  expect_equal(m$df, 6)
  expect_equal(m$df.residual, n - 6)
})

test_that("setup() sets class attribute", {
  expect_s3_class(m, "lmls")
})

# Tests if chol_info is created correctly -------------------------------------
test_that("chol_info created correctly:", {
  m <- init_beta(m)
  m <- init_gamma(m)
  m$chol_info_beta <- chol(info_beta(m))
  m$chol_info_gamma <- chol(info_gamma(m))

  expect_true(is.numeric(m$chol_info_gamma))
  expect_equal(dim(m$chol_info_gamma), c(3, 3))

  expect_true(is.numeric(m$chol_info_gamma))
  expect_equal(dim(m$chol_info_gamma), c(3, 3))
})

# Test data argument of lasso_lmls() ------------------------------------------
test_that("lasso_lmls() argument 'data' works with data frame", {
  dat <- data.frame(y = y, x4 = x1, x5 = x2, x6 = x3)
  m <- lasso_lmls(y ~ x4 + x6,
             ~ x5 + x6,
             data = dat,
             const = 1e-06)

  expect_error(m, NA)
})

test_that("lasso_lmls() argument 'data' works with list", {
  dat <- list(y = y, x4 = x1, x5 = x2, x6 = x3)
  m <- lasso_lmls(y ~ x4 + x6,
             ~ x5 + x6,
             data = dat,
             const = 1e-06)

  expect_error(m, NA)
})

# Testing init and update beta/gamma:
set.seed(969898)

n <- 1000
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
y <- rnorm(n, 0 + 1 * x1 + 1 * x3, exp(-3 + 1 * x2 + 1 * x3))
m <- lasso_lmls(y ~ x1 + x3,
                ~ x2 + x3,
                data = c(y, x1, x2, x3))

# Setting coefficients --------------------------------------------------------
beta <- gamma <- array(c(0, 0, 0), dim = c(3, 1))

m <- set_coef(m, "location", beta)
m <- set_coef(m, "scale", gamma)

test_that("setting coefficients works", {
  dimnames(beta) <- list(c("(Intercept)", "x1", "x3"))
  dimnames(gamma) <- list(c("(Intercept)", "x2", "x3"))

  expect_equal(coef(m, "location"), beta)
  expect_equal(coef(m, "scale"), gamma)
})

test_that("fitted values are updated", {
  expect_exactly(fitted(m, "location"), rep(0, n))
  expect_exactly(fitted(m, "scale"), rep(1, n))
})

test_that("residuals are updated", {
  expect_exactly(residuals(m), y)
})

# Coefficient matrices --------------------------------------------------------
test_that("coefficient matrix works", {
  cm <- coefmat(m, "location")

  expect_true(is.numeric(cm))
  expect_equal(dim(cm), c(3, 4))

  cn <- colnames(cm)

  expect_equal(rownames(cm), c("(Intercept)", "x1", "x3"))
  expect_equal(cn, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
})
