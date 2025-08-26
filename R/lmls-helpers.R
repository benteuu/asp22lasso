# log-likelihood --------------------------------------------------------------

#' @importFrom stats dnorm fitted

loglik <- function(m) {
  if (m$lasso){
    y <- m$y
    location <- fitted(m, "location")
    scale <- fitted(m, "scale")
    beta <- coef(m, "location")
    gamma <- coef(m, "scale")
    K_x <- if (ncol(m$x)>1) diag(c(0, rep(1, ncol(m$x)-1))) else matrix(0, 1, 1)
    K_z <- if (ncol(m$z)>1) diag(c(0, rep(1, ncol(m$z)-1))) else matrix(0, 1, 1)

    drop(sum(dnorm(y,
              location ,
              scale ,
              log = TRUE)) -
      m$lambda$location * sqrt(t(beta) %*% K_x %*% beta + m$const) -
      m$lambda$scale * sqrt(t(gamma) %*% K_z %*% gamma + m$const))
  } else {
    y <- m$y
    location <- fitted(m, "location")
    scale <- fitted(m, "scale")

    sum(dnorm(y, location, scale, log = TRUE))
  }
}

neg_loglik <- function(m, pred, y_test) {
  if (m$lasso){
    y <- y_test
    location <- pred$location
    scale <- pred$scale
    beta <- coef(m, "location")
    gamma <- coef(m, "scale")
    K_x <- if (ncol(m$x)>1) diag(c(0, rep(1, ncol(m$x)-1))) else matrix(0, 1, 1)
    K_z <- if (ncol(m$z)>1) diag(c(0, rep(1, ncol(m$z)-1))) else matrix(0, 1, 1)

    drop(-(sum(dnorm(y,
              mean = location,
              sd = scale,
              log = TRUE)) -
      m$lambda$location * sqrt(t(beta) %*% K_x %*% beta + m$const) -
      m$lambda$scale * sqrt(t(gamma) %*% K_z %*% gamma + m$const))
    )
  } else {
    y <- m$y
    location <- fitted(m, "location")
    scale <- fitted(m, "scale")

    - sum(dnorm(y, location, scale, log = TRUE))
  }
}

# derivatives -----------------------------------------------------------------

#' @importFrom stats fitted resid

score_beta <- function(m) {
  if (m$lasso){
    K_x <- if (ncol(m$x)>1) diag(c(0, rep(1, ncol(m$x)-1))) else matrix(0, 1, 1)
    beta <- coef(m, "location")

    drop(drop(resid(m, "response") / fitted(m, "scale")^2) %*% m$x) -
      (m$lambda$location * drop((t(beta) %*% K_x %*% beta + m$const)^(-1/2)) * drop(t(beta) %*% K_x))
  } else {
    drop((resid(m, "response") / fitted(m, "scale")^2) %*% m$x)
  }
}

#' @importFrom stats resid

score_gamma <- function(m) {
  if (m$lasso){
    K_z <- if (ncol(m$z)>1) diag(c(0, rep(1, ncol(m$z)-1))) else matrix(0, 1, 1)
    gamma <- coef(m, "scale")

    drop(drop(resid(m, "deviance")^2 - 1) %*% m$z) -
      (m$lambda$scale * drop((t(gamma) %*% K_z %*% gamma + m$const)^(-1/2)) * drop(t(gamma) %*% K_z))
  } else {
    drop((resid(m, "deviance")^2 - 1) %*% m$z)
  }
}

#' @importFrom stats fitted

info_beta <- function(m) {
  if (m$lasso){
    K_x <- if (ncol(m$x)>1) diag(c(0, rep(1, ncol(m$x)-1))) else matrix(0, 1, 1)
    beta <- coef(m, "location")

    crossprod(m$x/ drop(fitted(m, "scale"))) +
      m$lambda$location * (
        drop((t(beta) %*% K_x %*% beta + m$const)^(-0.5)) * K_x -
        drop((t(beta) %*% K_x %*% beta + m$const)^(-1.5)) * (K_x %*% beta) %*% t(K_x %*% beta)
      )
  } else {
    crossprod(m$x / fitted(m, "scale"))
  }
}

info_gamma <- function(m) {
  if (m$lasso){
    K_z <- if (ncol(m$z)>1) diag(c(0, rep(1, ncol(m$z)-1))) else matrix(0, 1, 1)
    gamma <- coef(m, "scale")

    2 * crossprod(m$z) +
      m$lambda$scale * (
        drop((t(gamma) %*% K_z %*% gamma + m$const)^(-0.5)) * K_z -
        drop((t(gamma) %*% K_z %*% gamma + m$const)^(-1.5)) * (K_z %*% gamma) %*% t(K_z %*% gamma)
      )
  } else {
    2 * crossprod(m$z)
  }
}

# helpers ---------------------------------------------------------------------

set_beta <- function(m, beta) {
  m$coefficients$location[] <- beta
  m$fitted.values$location <- drop(m$x %*% beta)
  m$residuals <- m$y - fitted(m, "location")
  m
}

set_gamma <- function(m, gamma) {
  m$coefficients$scale[] <- gamma
  m$fitted.values$scale <- exp(drop(m$z %*% gamma))
  m
}

set_coef_funs <- list(
  location = set_beta,
  scale = set_gamma
)

set_coef <- function(m, predictor, coef) {
  set_coef_funs[[predictor]](m, coef)
}

score_funs <- list(
  location = score_beta,
  scale = score_gamma
)

score <- function(m, predictor) {
  score_funs[[predictor]](m)
}

chol_info_funs <- list(
  location = function(m) chol(info_beta(m)),
  scale = function(m) m$chol_info_gamma
)

chol_info <- function(m, predictor) {
  chol_info_funs[[predictor]](m)
}

# multivariate normal distribution --------------------------------------------

#' @importFrom stats rnorm

rmvnorm <- function(n, mu = 0, chol_sig_inv) {
  dim <- nrow(chol_sig_inv)

  std_norm <- matrix(rnorm(dim * n), dim, n)
  scaled <- backsolve(chol_sig_inv, std_norm)
  shifted <- scaled + mu

  shifted
}

#' @importFrom stats dnorm

dmvnorm <- function(x, mu = 0, chol_sig_inv, log = FALSE) {
  std_norm <- drop(chol_sig_inv %*% (x - mu))
  correction <- sum(log(diag(chol_sig_inv)))

  log_prob <- dnorm(std_norm, log = TRUE)

  if (is.matrix(log_prob)) {
    log_prob <- colSums(log_prob) + correction
  } else {
    log_prob <- sum(log_prob) + correction
  }

  if (log) {
    log_prob
  } else {
    exp(log_prob)
  }
}

# summary helpers -------------------------------------------------------------

#' @importFrom generics tidy

coefmat <- function(m, predictor) {
  m <- tidy(m, predictor)
  out <- as.matrix(m[c("estimate", "std.error", "statistic", "p.value")])
  colnames(out) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  rownames(out) <- m$term
  out
}

#' @importFrom stats quantile

coefmat_samples <- function(m, predictor, type) {
  samples <- m[[type]][[predictor]]

  coefmat <- apply(samples, 2, function(x) {
    c(Mean = mean(x), quantile(x, c(0.025, 0.5, 0.975)))
  })

  t(coefmat)
}
