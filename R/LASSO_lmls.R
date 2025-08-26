#' @importFrom stats coef fitted lm.fit resid

init_beta <- function(m) {
  if (m$lasso) {
    K <- diag(c(0, rep(1, ncol(m$x)-1)))

    beta <- solve(t(m$x) %*% m$x + m$lambda$location * K) %*% t(m$x) %*% m$y

    m$coefficients$location <- beta
    m$fitted.values$location <- m$x %*% beta
    m$residuals <- m$y - m$fitted.values$location
  } else {
    fit <- lm.fit(m$x, m$y)

    m$coefficients$location <- coef(fit)
    m$fitted.values$location <- fitted(fit)
    m$residuals <- resid(fit)
  }
  m
}

#' @importFrom stats coef fitted lm.fit resid

init_gamma <- function(m) {
  # -(log(2) + digamma(1 / 2)) / 2 = 0.6351814
  # see https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
  # and https://arxiv.org/abs/1503.06266

  if (m$lasso) {
    s <- 0.6351814 + log(abs(resid(m, "response")))
    K <- if (ncol(m$z)>1) diag(c(0, rep(1, ncol(m$z)-1))) else matrix(0, 1, 1)

    gamma <- solve(t(m$z) %*% m$z + m$lambda$scale * K, t(m$z) %*% s)

    m$coefficients$scale <- gamma
    m$fitted.values$scale <- exp(m$z %*% gamma)
  } else {
    fit <- lm.fit(m$z, 0.6351814 + log(abs(resid(m, "response"))))

    m$coefficients$scale <- coef(fit)
    m$fitted.values$scale <- exp(fitted(fit))
  }
  m
}

#' @importFrom stats coef fitted lm.wfit resid

update_beta <- function(m) {
  if (m$lasso) {
    step <- backsolve(
      r = m$chol_info_beta,
      x = forwardsolve(
        l = m$chol_info_beta,
        x = drop(score_beta(m)),
        upper.tri = TRUE,
        transpose = TRUE
      )
    )
    m <- set_beta(m, coef(m, "location") + step)
  } else {
    fit <- lm.wfit(m$x, m$y, fitted(m, "scale")^(-2))

    m$coefficients$location <- coef(fit)
    m$fitted.values$location <- fitted(fit)
    m$residuals <- resid(fit)
  }
  m
}

#' @importFrom stats coef

update_gamma <- function(m) {
  if (m$lasso){
    step <- backsolve(
      r = m$chol_info_gamma,
      x = forwardsolve(
        l = m$chol_info_gamma,
        x = drop(score_gamma(m)),
        upper.tri = TRUE,
        transpose = TRUE
      )
    )
    m <- set_gamma(m, coef(m, "scale") + step)
  } else {
    step <- backsolve(
      r = m$chol_info_gamma,
      x = forwardsolve(
        l = m$chol_info_gamma,
        x = score_gamma(m),
        upper.tri = TRUE,
        transpose = TRUE
      )
    )
    m <- set_gamma(m, coef(m, "scale") + step)
  }
  m
}

#' @importFrom stats model.matrix update

setup <- function(location,
                  scale = ~1,
                  data = environment(location),
                  lasso = FALSE,
                  lambda = c(1, 1),
                  const,
                  light = TRUE,
                  call = NULL) {
  scale <- update(scale, paste(location[[2]], "~."))
  y <- eval(location[[2]], data, environment(location))
  x <- model.matrix(location, data)
  z <- model.matrix(scale, data)

  nobs <- length(y)
  df <- ncol(x) + ncol(z)

  m <- structure(
    list(
      y               = y,
      x               = x,
      z               = z,
      lasso           = lasso,
      lambda          = list(location = lambda[1], scale = lambda[2]),
      const           = const,
      nobs            = nobs,
      df              = df,
      df.residual     = nobs - df,
      light           = light,
      call            = call,
      term            = list(location = location, scale = scale),
      coefficients    = list(location = NULL, scale = NULL),
      fitted.values   = list(location = NULL, scale = NULL),
      residuals       = NULL,
      vcov            = list(location = NULL, scale = NULL),
      chol_info_beta  = NULL,
      chol_info_gamma = NULL,
      iterations      = NULL,
      loss            = NULL
    ),
    class = "lmls"
  )
  m
}

estimate <- function(m, maxit = 100, reltol = sqrt(.Machine$double.eps)) {
  m <- init_beta(m)
  m <- init_gamma(m)
  m$chol_info_beta <- chol(info_beta(m))
  m <- update_beta(m)

  it <- 0
  enough <- TRUE

  while (it < maxit && enough) {
    it <- it + 1

    before <- loglik(m)
    m$chol_info_gamma <- chol(info_gamma(m))
    m$chol_info_beta <- chol(info_beta(m))
    m <- update_gamma(m)
    m <- update_beta(m)
    after <- loglik(m)

    enough <- abs(after - before) > reltol * (abs(before) + reltol)
  }

  if (enough) {
    warning("Estimation did not converge, maximum number of iterations reached")
  }
  m$iterations <- it
  m
}

finish <- function(m) {
  m$vcov$location <- chol2inv(chol(info_beta(m)))
  m$vcov$scale <- chol2inv(chol(info_gamma(m)))

  if (m$light) {
    m$x <- m$z <- m$chol_info_gamma <- NULL
  }
  m
}

#' Gaussian location-scale regression
#'
#' @description
#'
#' The location-scale regression model assumes a normally distributed response
#' variable with one linear predictor for the mean (= the location) and one for
#' the standard deviation (= the scale). The standard deviation is mapped to
#' the linear predictor through a log link.
#'
#' This function sets up the model object and estimates it with maximum
#' likelihood.
#'
#' @param location A two-sided formula with the response variable on the LHS
#'                 and the predictor for the mean on the RHS.
#' @param scale A one-sided formula with the predictor for the standard
#'              deviation on the RHS.
#' @param data A data frame (or list or environment) in which to evaluate
#'             the `location` and `scale` formulas.
#' @param lasso A logical values that implies weather the model will be penalized with
#'              a L1 regularization.
#' @param lambda A list of values that are used as the lambdas in the penalization
#'               for location and scale coefficients.
#' @param const A small constant used in the quadratic approximation for the
#'              penalization.
#' @param light If `TRUE`, the design matrices are removed from the estimated
#'              model to save some memory.
#' @param maxit The maximum number of iterations of the Fisher scoring
#'              algorithm.
#' @param reltol The relative convergence tolerance of the Fisher scoring
#'               algorithm.
#'
#' @return
#'
#' A fitted linear model for location and scale as an `lmls` S3 object.
#' The object has at least the following entries:
#'
#' - `y`: the response vector
#' - `nobs`: the number of observations
#' - `df`: the degrees of freedom
#' - `df.residual`: the residual degrees of freedom
#' - `coefficients`: the regression coefficients as a list with the names
#'   `location` and `scale`
#' - `fitted.values`: the fitted values as a list with the names `location`
#'   and `scale`
#' - `residuals`: the response residuals
#' - `coefficients`: the variance-covariance matrices of the regression
#'   coefficients as a list with the names `location` and `scale`
#' - `iterations`: the number of iterations the Fisher scoring algorithm
#'   took to converge
#'
#' @examples
#' library(asp22lasso)
#' n <- 1000
#' x1 <- runif(n)
#' x2 <- runif(n)
#' x3 <- runif(n)
#' y <- rnorm(n, 0 + 1 * x1 + 1 * x3, exp(-3 + 1 * x2 + 1 * x3))
#' data <- data.frame(matrix(c(y, x1, x2, x3), ncol = 4))
#' names(data) <- c("y", "X1", "X2", "X3")
#' m <- lmls(y ~ X1 + X3, ~ X2 + X3, data = data)
#' summary(m)
#' plot(m)
#' qqnorm(m)
#' @export

lmls <- function(location,
                 scale = ~1,
                 data = environment(location),
                 lasso = FALSE,
                 lambda = NULL,
                 const = NULL,
                 light = TRUE,
                 maxit = 100,
                 reltol = sqrt(.Machine$double.eps)) {
  m <- setup(location, scale, data, lasso, lambda, const = const, light, match.call())
  m <- estimate(m, maxit, reltol)
  m <- finish(m)
  m
}

#' Cross validation function for a LASSO model
#'
#' @description
#'
#' The cross validation function uses the location-scale regression model and
#' evaluates for k-folds and a set of penalization terms lambda
#' (for each location and scale) the loss (negative log likelihood) for each set
#' of parameters. The combination of lambdas that over all folds yields the
#' smallest aggregated loss will be chosen as the optimal lambda parameters.
#' These are used to fit the returned model.
#'
#' @param location A two-sided formula with the response variable on the LHS
#'                 and the predictor for the mean on the RHS.
#' @param scale A one-sided formula with the predictor for the standard
#'              deviation on the RHS.
#' @param data A data frame (or list or environment) in which to evaluate
#'             the `location` and `scale` formulas.
#' @param k_fold The amount of folds data is split into. Each fold is once used
#'               as validation set while the rest serves as training data. A higher
#'               number of folds increases the computational time but yields more
#'               accurate results for the loss.
#' @param steps The number of lambda values that the model will be run on. A higher
#'              number of steps will increase the computational time but yields more
#'              accurate results for the optimal lambda values.
#' @param const Refers to the constant in lmls() that is used in the quadratic
#'              approximation for the penalty term.
#'
#' @import broom dplyr glmnet pracma readr
#'
#' @return
#'
#' A fitted linear model for location and scale as an `lmls` S3 object.
#' The object has at least the following entries:
#'
#' - `y`: the response vector
#' - `x`: design matrix for the location variables
#' - `z`: design matrix for the scale variables
#' - `lasso`: binary, TRUE when model is lasso
#' - `lambda`: optimal lambdas as a list with the names `location`
#'   and `scale`
#' - `const`: small constant that is used to approximate the
#'   penalty term
#' - `nobs`: the number of observations
#' - `df`: the degrees of freedom
#' - `df.residual`: the residual degrees of freedom
#' - `call`: call of the model
#' - `term`: model specifications as a list with the names `location`
#'   and `scale`
#' - `coefficients`: the regression coefficients as a list with the names
#'   `location` and `scale`
#' - `fitted.values`: the fitted values as a list with the names `location`
#'   and `scale`
#' - `residuals`: the response residuals
#' - `vcov`: the variance-covariance matrices of the regression
#'   coefficients as a list with the names `location` and `scale`
#' - `chol_info_gamma`: cholesky decomposition for info_gamma(m)
#' - `chol_info_beta`: cholesky decomposition for info_beta(m)
#' - `iterations`: the number of iterations the Fisher scoring algorithm
#'   took to converge
#' - `loss`: the minimal loss that corresponds to optimal lambdas
#' - `loss matrix`: Outputs the loss values for each combination of lambdas
#'    as a matrix. The rows are the location lambdas and the columns are the
#'    scale lambdas.
#'
#' @examples
#' library(asp22lasso)
#' n <- 1000
#' x1 <- runif(n)
#' x2 <- runif(n)
#' x3 <- runif(n)
#' y <- rnorm(n, 0 + 1 * x1 + 1 * x3, exp(-3 + 1 * x2 + 1 * x3))
#' data <- data.frame(matrix(c(y, x1, x2, x3), ncol = 4))
#' names(data) <- c("y", "X1", "X2", "X3")
#' m <- lasso_lmls(y ~ X1 + X3, ~ X2 + X3, data = data)
#' summary(m)
#' plot(m)
#' qqnorm(m)
#' @export

lasso_lmls <- function(location,
               scale = ~1,
               data,
               k_fold = 10,
               steps = 20,
               const = 1e-06
               ){
  call <- match.call()

  lambda_loc <- logspace(1.5, 2, n = steps)
  lambda_scale <- logspace(1.5, 2, n = steps)
  pred_error <- matrix(0, nrow = length(lambda_loc), ncol = length(lambda_scale))

  # 1) Randomly split the data set into k-subsets (or k-fold)
  data <- as.data.frame(data)
  datalist <- split(data, sample(1:k_fold, nrow(data), replace=T)) #list of k same-sized elements that are slices of the data
  loss <- matrix(0, nrow=length(lambda_loc), ncol=length(lambda_scale))

  # 2) Reserve one subset and train the model on all other subsets
  for (i in 1:k_fold) {
    #split data in k folds
    data_val <- datalist[[i]]    #ith of the k folds, validation set
    data_train <- datalist[-i]   #rest of the data without ith of the k folds, training set
    data_train <- bind_rows(data_train) #convert list to dataframe

    for (j in 1:length(lambda_loc)){
      for (l in 1:length(lambda_scale)){
        m <- lmls(location = location,
                  scale = scale,
                  data = data_train,
                  light = FALSE,
                  lasso = TRUE,
                  lambda = c(lambda_loc[j], lambda_scale[l]),
                  const = const
                  )
        pred <- predict.lmls(object = m, newdata = data_val, type = "response")
        y_test <- eval(location[[2]], data_val, environment(location))
        loss[j,l] <- loss[j,l] + neg_loglik(m, pred = pred, y_test = y_test)
      }
    }
  }
  # finding the best lambdas
  minimum  <- which(loss==min(loss), arr.ind = TRUE)
  opt_loc <- lambda_loc[minimum[1]]
  opt_scale <- lambda_scale[minimum[2]]

  m_opt <- lmls(location,
                scale,
                data,
                light = FALSE,
                lasso = TRUE,
                lambda = c(opt_loc, opt_scale),
                const = const)
  m_opt$loss <- loss[minimum[1], minimum[2]]
  m_opt$lossmatrix <- loss

  return(m_opt)
}
