########################################################
# 1.
# Write a function for the logistic regression.
########################################################

glm_function <- function(Y, X, data, stval){
  
  # (1) Log-Likelihood function
  logit.ll <- function(X,Y,theta){
    mu <- theta[1] + X %*% theta[2:length(theta)]
    
    l.link <- function(x){
      exp(x)/(1+exp(x))
    }
    logl <- sum(Y*log(l.link(mu)) + (1 - Y)*log(1-l.link(mu)))
    return(-logl)
  }
  
  # (2) define n and k
  n = nrow(data) # sample size
  k = ncol(X) # number of betas
  
  # (3) optimization
  res <- optim(par=stval, fn=logit.ll, Y=Y, X=X, hessian=TRUE)
  
  # (4) calculate the standard errors
  se <- sqrt(diag(solve(res$hessian))) # gives us the standarddeviation
  
  # (5) calculate z value
  z = res$par/se
  
  # (6) calculate p value
  p = 2*pt(abs(z), df=n-(k+1), lower.tail= FALSE)
  
  # (6a) add stars
  stars = ifelse(p<0.001,"***",
                 ifelse(p<0.01,"**", 
                        ifelse(p<0.05,"*",
                               ifelse(p<0.1,"."," "))))
  p = ifelse(p < 2*10^(-16), 
             c("< 2e-16"), format(signif(p, digits = 3)))
  
  # (7) generate a data frame
  result = as.data.frame(cbind(c("(Intercept)", colnames(data[2:(k+1)])),
                               sprintf("%.6f", round(res$par, digits=6)),
                               sprintf("%.6f", round(se, digits=6)),
                               sprintf("%.3f", round(z, digits=3)),
                               p))
  
  colnames(result) = c("Coefficients","Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  result$' ' = stars
  
  # (8) print the output
  print(result)
  cat(paste0("---",
             "\n",
             "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",
             "\n",
             "---",
             "\n"))
}
