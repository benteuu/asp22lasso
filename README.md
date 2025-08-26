
<!-- README.md is generated from README.Rmd. Please edit that file -->

# asp22lasso

<!-- badges: start -->
<!-- badges: end -->

`asp22lasso` is a R package which extends the existing [`lmls`
package](https://github.com/hriebl/lmls) with LASSO (Least Absolute
Shrinkage and Selection Operator) regression. `lmls` is a
multi-predictor regression with explanatory variables for the mean (=
the location) and the standard deviation (= the scale) of a normally
distributed response variable. LASSO regression adds a penalty to these
coefficients and shrinks them towards or to zero. The optimal
hyperparameter lambda, which controls the strength of the LASSO penalty,
is found via Cross-Validation.

## Installation

You can download the development version of asp22lasso by cloning the
[gitlab repository](https://gitlab.gwdg.de/asp22/asp22lasso) and then
install the package with *devtools::install( )*. The function requires
the file path to the downloaded package.

## Example

This is a basic example which shows you how to solve a common problem:

This code simulates some basic data.

``` r
library(asp22lasso)

set.seed(1346)

n <- 1000
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
y <- rnorm(n, 0 + 1 * x1 + 1 * x3, exp(-3 + 1 * x2 + 1 * x3))
```

To estimate a LASSO-location-scale linear model, the function
*lasso_lmls( )* is used. The input arguments location, scale and data
are set with the generated data. Additional parameters are not modified
for this example. It uses the default values (see the function’s
description for further input options). The function outputs the model
refit on the full dataset with the optimal hyperparameters found by the
cross-validation. <br> To view the regression output, the function
*summary( )* is used.

``` r
m_opt <- lasso_lmls(location = y ~ x1 + x3,
                scale = ~ x2 + x3,
                data = c(y, x1, x2, x3))

summary(m_opt)
#> 
#> Call:
#> lmls(location = location, scale = scale, data = data, lasso = TRUE, 
#>     lambda = c(opt_loc, opt_scale), const = const, light = FALSE)
#> 
#> Deviance residuals:
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> -3.285000 -0.668400  0.002406 -0.002902  0.653500  2.901000 
#> 
#> Chosen Lambda values (location/scale):
#>  12.58925 12.58925
#> 
#> Location coefficients (identity link):
#>              Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) -0.005254   0.008286  -0.634    0.526    
#> x1           1.006058   0.012533  80.270   <2e-16 ***
#> x3           1.007986   0.013834  72.865   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Scale coefficients (log link):
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) -2.93087    0.05696  -51.45   <2e-16 ***
#> x2           0.92247    0.07460   12.37   <2e-16 ***
#> x3           0.90465    0.07707   11.74   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual degrees of freedom: 994
#> Log-likelihood: 582.136
#> AIC: -1152.27
#> BIC: -1122.82
```
