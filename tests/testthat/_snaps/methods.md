# print() works

    Code
      print(m)
    Output
      
      Call:
      lmls(location = location, scale = scale, data = data, lasso = TRUE, 
          lambda = c(opt_loc, opt_scale), const = const, light = FALSE)
      
      Location coefficients (identity link):
                   [,1]     
      (Intercept)  -0.004321
      x1            1.004162
      x3            0.995498
      
      Scale coefficients (log link):
                   [,1]   
      (Intercept)  -2.8938
      x2            0.9337
      x3            0.8283
      

# summary() works

    Code
      summary(m)
    Output
      
      Call:
      lmls(location = location, scale = scale, data = data, lasso = TRUE, 
          lambda = c(opt_loc, opt_scale), const = const, light = FALSE)
      
      Deviance residuals:
           Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
      -3.826000 -0.645400 -0.016040 -0.004047  0.659600  3.740000 
      
      Chosen Lambda values:
      Location: 31.623 / Scale: 31.623
      
      Location coefficients (identity link):
                   Estimate Std. Error t value Pr(>|t|)    
      (Intercept) -0.004321   0.009245  -0.467     0.64    
      x1           1.004162   0.013008  77.195   <2e-16 ***
      x3           0.995498   0.013511  73.679   <2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      Scale coefficients (log link):
                  Estimate Std. Error t value Pr(>|t|)    
      (Intercept) -2.89381    0.06222  -46.51   <2e-16 ***
      x2           0.93370    0.07643   12.22   <2e-16 ***
      x3           0.82827    0.07467   11.09   <2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      Residual degrees of freedom: 994
      Log-likelihood: 502.119
      AIC: -992.238
      BIC: -962.791
      

