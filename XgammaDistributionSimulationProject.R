# Simulation Studies for Xgamma distribution


#_____________________________________________________________________________

# Loading Packages

library(MASS) # fitdistr function
library(tidyverse)

#_____________________________________________________________________________

### Function to generate Xgamma observations

"rxgamma" = function(n, rate)
{
  threshold = rate / (1+rate)
  UNIF = runif(n = n, min = 0, max = 1)
  EXP = rexp(n = n, rate = rate)
  GAMMA = rgamma(n = n, shape = 3, rate = rate)
  XGAMMA = ifelse(UNIF <= threshold, EXP, GAMMA) 
  return(XGAMMA)
}

#_____________________________________________________________________________

# Density function of Xgamma - x > 0, rate > 0

"dxgamma" = function(x, rate)
{
  (rate^2 / (1 + rate))*(1 + (rate/2)*(x^2))*(exp(-(rate*x)))
}

#_____________________________________________________________________________

# Plot density function of Xgamma for selected values of theta

# Define the range of x values
X = seq(0, 10, length.out = 500)

# Specify rate values
THETAS.fx = c(0.5, 1, 5, 10) 

# Density values for all rates using dxgamma function
Y = sapply(THETAS.fx, 
           function(rate) dxgamma(X, rate=rate))

# Adjust plot margins
par(mar=c(4.5, 4, 1, 1))

# Plot the density functions
matplot(X, Y, 
        type = "l", col = rainbow(length(THETAS.fx)), 
        lty = 1:4, lwd = 2,
        xlab = "x", ylab = "f(x)", ylim = c(0, 0.5), 
        main = "")

# Add a legend
legend("topright", 
       legend = paste("theta =", THETAS.fx), 
       col = rainbow(length(THETAS.fx)), 
       lty = 1:4, lwd = 2, bty = "n")


#_____________________________________________________________________________

# Log-Likelihood Function for Xgamma distribution

"loglik.xgamma" = function(rate, x)
{
  n = length(x)
  loglik = 2*n*log(rate) - n*log(1 + rate) + 
    sum(log(1 + (rate/2)*(x^2))) - rate*sum(x)
  
  return(-loglik) 
}

#_____________________________________________________________________________

# Function to compute Maximum Likelihood Estimate(MLE) theta.hat of theta

"xgamma.theta.ML" = function(x)
{
  theta.optim = optim(par=0.1, 
                      fn = loglik.xgamma,
                      x = x,
                      method = "Brent", 
                      lower = 0, upper = 5000,
                      hessian = FALSE)
  return(theta.optim$par)
}

"exp.theta.ML" = function(x)
{
  mle = fitdistr(
    x = x, 
    densfun = "exponential"
  )
  return(mle$estimate)
}

#_____________________________________________________________________________

### MSE and Bias function (Secondary function)

"Bias.MSE" = function(est, param)
{
  # est = vector of estimates
  # param = parameter to estimate
  Bias = mean(est - param)
  MSE = mean((est - param)^2)
  tab = round(c(Bias, MSE), 5)
  names(tab) = c("Bias", "MSE")
  return(tab)
}


#_____________________________________________________________________________
  
# Monte Carlo Simulation function
# Computes Bias and MSE 

"sim.xgamma.exp" = function(N, n, rate)
{
  
  #________________________________________________________________
  
  # Generate random samples from Xgamma distribution 
  samples = replicate(
    n = N,
    rxgamma(n=n, rate = rate)
  )
  
  #________________________________________________________________
  
  # Maximum Likelihood Estimates for samples
  xgamma.theta.hat = apply(
    samples,
    MARGIN = 2,
    function(X) xgamma.theta.ML(X)
  )
  
  # Bias and MSE for xgamma theta.hat
  xgamma.res = Bias.MSE(est = xgamma.theta.hat, param = rate)

  #________________________________________________________________

  
  # Maximum Likelihood Estimates for exponential samples
  exp.theta.hat = apply(
    samples,
    MARGIN = 2,
    function(X) exp.theta.ML(X)
  )
  
  # Compute Bias and MSE
  exp.res = Bias.MSE(est = exp.theta.hat, param = rate)
  
  #________________________________________________________________
  
  tab = rbind(xgamma.res, exp.res)
  row.names(tab) = c("Xgamma", "Exponential")
  
  return(tab)
}


#_____________________________________________________________________________

# Testing Simulation function

R = 10000
n = 20
rate = 0.1
# sim.xgamma.exp(N=N, n=n, rate = rate)

#_____________________________________________________________________________

# Sample sizes and parameter values objects

N = c(20, 40, 100) # Sample sizes

THETA = c(0.1, 0.5, 
      1.0, 1.5, 3, 6) # Parameter values

# Generate combinations of sample sizes and parameters

size_param_combos = expand.grid(n = N, rate = THETA)

size_param_combos

#_____________________________________________________________________________

# Presenting simulation results

sim.res = mapply(sim.xgamma.exp, R, size_param_combos$n, size_param_combos$rate, SIMPLIFY = "array")
dimnames(sim.res)[[3]] = paste("n_", size_param_combos$n, "Theta_", size_param_combos$rate, sep="")
sim.res
