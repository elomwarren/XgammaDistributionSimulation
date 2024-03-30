# Simulation Studies for Xgamma distribution

# library(tidyr) # presentation of data
# library(dplyr) # data manipulation
# library(MASS) # for fitdistr
pacman::p_load(tidyr, dplyr, MASS)



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

# Sample sizes and parameter values objects

N = c(20, 40, 100) # Sample sizes

THETA = c(0.1, 0.5, 
      1.0, 1.5, 3, 6) # Parameter values

# Generate combinations of sample sizes and parameters
size_param_combos = expand.grid(n = N, rate = THETA)

size_param_combos

### Generate Xgamma observations for all combinations of sample size and parameter rate


list_Data = lapply(
  1:nrow(size_param_combos), 
  function(i) # Creating an anonymous function 
  {
    rxgamma(n = size_param_combos$n[i], 
            rate = size_param_combos$rate[i])
  }
  )

# Create names based on combinations of n and rate
combination_names = paste("n", size_param_combos$n, 
                          "_rate", size_param_combos$rate, 
                          sep = "")

# Set names for the list elements
names(list_Data) = combination_names


# Density function of Xgamma
# x > 0, rate > 0
"dxgamma" = function(x, rate)
{
  (rate^2 / (1 + rate))*(1 + (rate/2)*(x^2))*(exp(-(rate*x)))
}


# Plot density function of Xgamma for selected values of theta

# Define the range of x values
X = seq(0, 10, length.out = 500)

# Specify rate value
THETAS.fx = c(0.5, 1, 5, 10) 

# Density values for all rates using dxgamma function
Y = sapply(THETAS.fx, 
           function(rate) dxgamma(X, rate=rate))

# Plot the density functions
matplot(X, Y, 
        type = "l", col = rainbow(length(THETAS.fx)), 
        lty = 1:4, lwd = 2,
        xlab = "x", ylab = "f(x)", ylim = c(0, 0.5), 
        main = "Xgamma Density Functions")

# Add a legend
legend("topright", 
       legend = paste("theta =", THETAS.fx), 
       col = rainbow(length(THETAS.fx)), 
       lty = 1:4, lwd = 2, bty = "n")


## Log-Likelihood Function for Xgamma

"loglik.xgamma" = function(rate, x)
{
  n = length(x)
  loglik = 2*n*log(rate) - n*log(1 + rate) + 
    sum(log(1 + (rate/2)*(x^2))) - rate*sum(x)
  
  return(-loglik) 
}


## Function to compute Maximum Likelihood Estimate MLE theta.hat of theta
"theta.ML" = function(x)
{
  theta.optim = optim(par=0.1, 
                      fn = loglik.xgamma,
                      method = "Brent", 
                      lower = 0, upper = 5000,
                      hessian = FALSE)
  return(theta.optim$par)
}

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


## Monte Carlo Simulation function for Xgamma
# Generates Bias and MSE of the estimator theta.hat

"simple.Sim.xgamma" = function(N, n, rate)
{
  # Generate samples from Xgamma distribution 
  samples = sapply(
    1:N,
    function(i) rxgamma(n=n, rate = rate)
  )
  
  # Maximum Likelihood 
  res.rate.hat = apply(
    samples, 
    MARGIN = 2, # function applied to columns
    function(x)
    {
      optim(par = 0.2,
            fn = loglik.xgamma, 
            x = x, 
            method = c("Brent"), 
            lower = 0, upper = 3000,
            hessian = T
            )
    }
  )
  
  # Finding rate.hat
  rate.hat = sapply(
    1:N,
    function(x)res.rate.hat[[x]]$par
  )
  
  # Bias and MSE computation
  tab = Bias.MSE(rate.hat = rate.hat, rate = rate)

  return(tab)
}


# Testing Simulation function

N = 10000
n = 20
rate = 0.1
simple.Sim.xgamma(N=N, n=n, rate = rate)




# Simulation for n = 20
N = 10000
n = 20
Bias_MSE_mat = data.frame(
  t(sapply(THETA, 
           function(rate)
             simple.Sim.xgamma(N=N, n=n, rate = rate))))



### Calling simulation function for selected values of n and rate

n.rate = expand.grid(rate = THETA, n = N)
N = rep(10000, length(n.rate))



### Simulation function for all values of n and rate

"Sim.xgamma" = function(N, n, rate)
{
  Bias.MSE_matrix = data.frame(
    t(
      mapply(simple.Sim.xgamma, N=N, n=n, rate = rate)), 
    row.names = paste("n", n, "_rate", rate, sep = ""))

  return(Bias.MSE_matrix)
}


## Monte Carlo Simulation function for exponential

# Simulation function for exponential distribution

"simple.Sim.exp" = function(N, n, rate)
{
  # Generate random data from exponential distribution
  exp.Data = rexp(n = n, rate = rate)
  
  # Maximum Likelihood estimation
  res.rate.hat = fitdistr(
    x = exp.Data, 
    densfun = "exponential"
  )
  
  # Return MLE
  rate.hat = res.rate.hat$estimate
  
  # Compute Bias and MSE
  tab = Bias.MSE(rate.hat = rate.hat, rate = rate)
  
  return(tab)
}


# Testing exponential distribution simulation function

N = 10000
n = 20
rate = 0.1
simple.Sim.exp(N=N, n=n, rate = rate)


### Main simulation function for exponential

"Sim.exp" = function(N, n, rate)
{
  Bias.MSE_matrix = data.frame(
    t(
      mapply(simple.Sim.exp, N=N, n=n, rate = rate)), 
    row.names = paste("n", n, "_rate", rate, sep = ""))

  return(Bias.MSE_matrix)
}



## Simulations results - Bias and MSE 

# For Xgamma

Bias.MSE.xgamma = Sim.xgamma(N, n.rate$n, n.rate$rate)

Bias.MSE.xgamma.df = data.frame(
  n = rep(N, each = 6),
  rate = rep(THETA, times = 3),
  Bias = Bias.MSE.xgamma$Bias,
  MSE = Bias.MSE.xgamma$MSE
)

Sim.Bias.MSE.xgamma = Bias.MSE.xgamma.df %>%
  pivot_wider(names_from = n, 
              values_from = c(Bias, MSE),
              names_vary = "slowest")




# For Exponential

Bias.MSE.exp = Sim.exp(N, n.rate$n, n.rate$rate)

Bias.MSE.exp.df = data.frame(
  n = rep(N, each = 6),
  rate = rep(THETA, times = 3),
  Bias = Bias.MSE.exp$Bias,
  MSE = Bias.MSE.exp$MSE
)

Sim.Bias.MSE.exp = Bias.MSE.exp.df %>%
  pivot_wider(names_from = n, 
              values_from = c(Bias, MSE),
              names_vary = "slowest")





## Presenting Simulations results Xgamma and Exponential

knitr::kable(Sim.Bias.MSE.xgamma, 
             caption = "Average bias and MSE of the MLE for Xgamma")

knitr::kable(Sim.Bias.MSE.exp, 
             caption = "Average bias and 
             MSE of the MLE for Exponential")

