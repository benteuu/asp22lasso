# Function to simulate data:
sim_data <- function(n, diff_mat = FALSE, p=0, seed = 42){
  #p is the correlation inbetween the variables
  set.seed(seed)
  num_beta <- num_gamma <- 8

  true_beta <- c(0, 0.15, 0.3, 0.5, 0.7, 1, 1.2, 1.5)
  true_gamma <- c(0, 0.03, 0.05, 0.1, 0.15, 0.19, 0.28, 0.31)

  true_gamma <- true_gamma / sum(true_gamma)

  #mean values for the X values
  means <- runif(min = -0.3, max = 0.3, n = length(true_beta))

  #variances for the X values
  variances <- diag(runif(min = 0.01, max = 0.4, n = length(true_beta)), length(true_beta))

  #matrix for correlation p
  upper <- matrix(p, length(true_beta), length(true_beta))
  upper[lower.tri(upper, diag=TRUE)] <- 0
  lower <- matrix(p, length(true_beta), length(true_beta))
  lower[upper.tri(lower, diag=TRUE)] <- 0

  X <- matrix(mvtnorm::rmvnorm(n, mean = means, sigma = variances + lower + upper), nrow = n)

  if(diff_mat) {
    Z <- matrix(c(rep(0, n), rnorm(n = n*(length(true_gamma)-1), mean = 1, sd = 1)), nrow = n)
  } else {
    Z <- X
  }

  y <- rep(NA, n)

  for (i in 1:n){
    mean <- X[i,] %*% true_beta
    log_sd <- Z[i,] %*% true_gamma
    y[i] <- rnorm(n = 1, mean = mean, sd = exp(log_sd))
  }

  X <- data.frame(X)
  Z <- data.frame(Z)

  output <- list("true_beta" = true_beta, "true_gamma" = true_gamma,
                 "X" = X, "Z" = Z, "y" = y)
  return(output)
}

# Simulation Study:
#------------------------------------------------------------------------------
simulation <- function(nrun = 5, nobs = 1000, p = 0){
  require(dplyr)
  # bias of the coefficients
  bias <- data.frame(matrix(0, nrow = 8, ncol = 4))
  names(bias) <- c("lasso location", "lasso scale", "lmls location", "lmls scale")
  # coefficients
  coe <- array(0, dim = c(8, 4, nrun))
  # standard deviation
  s_d <- data.frame(matrix(0, nrow = 8, ncol = 4))
  names(s_d) <- c("lasso location", "lasso scale", "lmls location", "lmls scale")
  msd <- data.frame(matrix(0, nrow = 1, ncol = 4))
  names(msd) <- c("lasso location", "lasso scale", "lmls location", "lmls scale")
  # excluded variables percentage
  excl <- data.frame(matrix(0, nrow = 8, ncol = 4))
  # Lambdas
  lambdas <- data.frame(matrix(0, nrow = nrun, ncol = 2))
  names(lambdas) <- c("Lambda location", "Lambda scale")


  s <- 1345
  for (i in 1:nrun){
    print(c("Runde:", i))
    s <- s + 1

    data <- sim_data(nobs, seed = s, p = p)
    d <- cbind(data.frame(data$y), data.frame(data$X))
    names(d)[1] <- "y"

    # models
    m_opt <- lasso_lmls(location = y ~ X1 + X2 + X3 + X4+ X5 + X6 + X7 + X8,
                scale = ~ X1 + X2 + X3 + X4+ X5 + X6 + X7 + X8,
                data = d,
                k_fold = 10,
                steps = 20
                )
    m_opt_lmls <- lmls(location = y ~ X1 + X2 + X3 + X4+ X5 + X6 + X7 + X8,
                       scale = ~ X1 + X2 + X3 + X4+ X5 + X6 + X7 + X8,
                       data = d,
                       lasso = FALSE
                       )
    lambdas[i,1] <- m_opt$lambda$location
    lambdas[i,2] <- m_opt$lambda$scale
    # bias
    # for LASSO_lmls
    bias[,1] <- bias[,1] + (m_opt$coefficients$location[2:9] - data$true_beta)
    bias[,2] <- bias[,2] + (m_opt$coefficients$scale[2:9] - data$true_gamma)
    # for lmls
    bias[,3] <- bias[,3] + (m_opt_lmls$coefficients$location[2:9] - data$true_beta)
    bias[,4] <- bias[,4] + (m_opt_lmls$coefficients$scale[2:9] - data$true_gamma)

    # saving the coefficients
    for (j in 1:8){
      # j + 1 because we want to skip the intercept
      coe[j,1,i] <- m_opt$coefficients$location[j+1]
      coe[j,2,i] <- m_opt$coefficients$scale [j+1]
      coe[j,3,i] <- m_opt_lmls$coefficients$location[j+1]
      coe[j,4,i] <- m_opt_lmls$coefficients$scale [j+1]
    }
  }
  # mean bias
  bias <- bias / nrun


  for (j in 1:8){
    # standard deviation
    s_d[j,1] <- sqrt(mean((coe[j,1,] - mean(coe[,1,]))^2))
    s_d[j,2] <- sqrt(mean((coe[j,2,] - mean(coe[,2,]))^2))
    s_d[j,3] <- sqrt(mean((coe[j,3,] - mean(coe[,3,]))^2))
    s_d[j,4] <- sqrt(mean((coe[j,4,] - mean(coe[,4,]))^2))
  }

  for (k in 1:nrun){
    for (j in 1:8){
      # for every coefficient check if it's excluded
      # LASSO_lmls
      excl[j,1] <- excl[j,1] + between(0,
                                       coe[j,1,k] - s_d[j, 1],
                                       coe[j,1,k] + s_d[j, 1])
      excl[j,2] <- excl[j,2] + between(0,
                                       coe[j,2,k] - s_d[j, 2],
                                       coe[j,2,k] + s_d[j, 2])
      # for lmls
      excl[j,3] <- excl[j,3] + between(0,
                                       coe[j,3,k] - s_d[j, 3],
                                       coe[j,3,k] + s_d[j, 3])
      excl[j,4] <- excl[j,4] + between(0,
                                       coe[j,4,k] - s_d[j, 4],
                                       coe[j,4,k] + s_d[j, 4])
    }
  }

  # mean sd over all runs
  msd[1] <- mean(s_d[,1])
  msd[2] <- mean(s_d[,2])
  msd[3] <- mean(s_d[,3])
  msd[4] <- mean(s_d[,4])

  # exclusion rate: #of exclusion/ nrun
  excl <- excl / nrun

  # save true values
  true_vals <- list(data$true_beta, data$true_gamma)
  names(true_vals) <- c("True_Beta", "True_Gamma")

  results <- list(bias, s_d, msd, excl, lambdas, true_vals)
  names(results) <- c("Bias", "Standard Deviation", "Mean Standard Deviation",
                      "How often Excluded", "Lambas of Models", "True Coef")
  return(results)
}

#Saving the data:
library(here)
library(backports)

res1 <- simulation(nrun=100, nobs=1000, p=0)
saveRDS(res1, file = (paste0(here(), "/simulation/res1.RData")))
res2 <- simulation(nrun=100, nobs=5000, p=0)
saveRDS(res2, file = (paste0(here(), "/simulation/res2.RData")))


