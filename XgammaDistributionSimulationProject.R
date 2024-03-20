

# Simulation Studies for Xgamma distribution

# library(optimx) # advanced version of optim
# library(xtable)
# library(tidyr) # presentation of data
# library(dplyr) # data manipulation
# library(MASS) # fitdistr
pacman::p_load(tidyr, dplyr, MASS, optimx, xtable)



## Generating Observations from the Xgamma distribution

# Secondary function to choose which value to pick

# "CHOOSE" = function(u, x, y, thres)
# {
#   if (u <= thres)
#   {
#     return(x)
#   }
#   else
#   {
#     return(y)
#   }
# }



### Function to generate random XGamma observations

"rxgamma" = function(n, Theta)
{
  UNIF = runif(n = n, min = 0, max = 1)
  EXP = rexp(n = n, rate = Theta)
  GAMMA = rgamma(n = n, shape = 3, scale = 1/Theta)
  
    "CHOOSE" = function(u, x, y, thres)
  {
    if (u <= thres)
    {
      return(x)
    }
    else
    {
      return(y)
    }
  }
  
  Obs = mapply(CHOOSE, UNIF, EXP, GAMMA, thres = Theta/(1+Theta))
  return(Obs)
}




# Setting up sample sizes and parameter values

# Sample sizes
n_values = c(20, 40, 100)
# Parameter values
Theta_values = c(0.1, 0.5, 1.0, 1.5, 3, 6)

# Generate combinations of sample sizes and parameters
size_param_combos = expand.grid(n = n_values, Theta = Theta_values)

size_param_combos



### Generate Xgamma observations for all combinations of sample size and parameter Theta


list_Data = lapply(
  1:nrow(size_param_combos), 
  function(i) # Creating an anonymous function 
  {
    rxgamma(n = size_param_combos$n[i], 
            Theta = size_param_combos$Theta[i])
  }
  )

# Create names based on combinations of n and Theta
combination_names = paste("n", size_param_combos$n, 
                          "_Theta", size_param_combos$Theta, 
                          sep = "")

# Set names for the list elements
# setNames(list_Data, combination_names) # stopped working?!

names(list_Data) = combination_names




# HISTOGRAMS of Data

# HIST = lapply(1:length(list_Data),
#               function(i)
#               {
#                 hist(list_Data[[i]],
#                      xlab = paste("Observations of Data",
#                                   names(list_Data)[i]), main = "")
#               }
# )

# Histogram for dataset with n=20 and Theta = 0.1
# hist(list_Data[[1]], 
#      xlab = paste("Observations of Data", 
#                   names(list_Data)[1]), main = "")



# DENSITY Graph of Data

# PLOT DENSITY

# DENSITY = lapply(
#   1:length(list_Data),
#   function(i)
#   {
#     plot(density(list_Data[[i]]),
#          xlab = paste("Density" ,
#                       names(list_Data)[i]), main = "")
#   }
#   )

# Density for parameter combinations with size 20

DENSITY_n20 = lapply(
  c(1, 4, 7, 10, 13, 16),
  function(i)
  {
     plot(density(list_Data[[i]]),
          xlab = paste("Density" ,
                       names(list_Data)[i]), main = "")
  }
)




## Probability density function pdf of the Xgamma distribution


### Function for pdf of Xgamma

# x > 0, Theta > 0
"den_xgamma" = function(Theta, x)
{
  (Theta^2 / (1 + Theta))*(1 + (Theta/2)*(x^2))*(exp(-(Theta*x)))
}



# Graph pdf of Xgamma(Theta) for selected values of Theta

# Define the range of x values
x_values = seq(0, 10, length.out = 1000)

# Specify Theta value
Theta_values_pdf = c(0.5, 1, 5, 10) 

# Density values for all Thetas using den_xgamma function
density_values = sapply(Theta_values_pdf, 
                        function(Theta) den_xgamma(Theta, x_values))

# Plot the density functions
matplot(x_values, density_values, 
        type = "l", col = rainbow(length(Theta_values_pdf)), 
        lty = 1:4, lwd = 2,
        xlab = "x", ylab = "f(x)", ylim = c(0, 0.5), 
        main = "Xgamma Density Functions")

# Add a legend
legend("topright", 
       legend = paste("Theta =", Theta_values_pdf), 
       col = rainbow(length(Theta_values_pdf)), 
       lty = 1:4, lwd = 2, bty = "n")





## Log-Likelihood Function for Xgamma


"loglik_Xgamma" = function(Theta, x)
{
  n = length(x)
  loglik = 2*n*log(Theta) - n*log(1 + Theta) + 
    sum(log(1 + (Theta/2)*(x^2))) - Theta*sum(x)
  
  return(-loglik) 
}



## Finding Maximum Likelihood Estimator MLE Theta.Hat of Theta


# vec_Theta.Hat = lapply(
#   1:length(list_Data),
#   function(i)
#   {
#     tryCatch(unlist(
#          optimx(par = 0.2,
#            fn = loglik_Xgamma,
#            method = "nlminb",
#            x = list_Data[[i]]
#     )),
#     error = function(e)rep(NA))
#   }
# )


# Try and Error on optim function to find MLE

vec_Theta.Hat = sapply(
  1:length(list_Data),
  function(i)
  {
    optim(par = 0.2,
          fn = loglik_Xgamma, 
          x = list_Data[[i]], 
          method = c("Brent"), 
          lower = 0, upper = 3000,
          hessian = T
          )
  }
)



# vec_Theta.Hat
# # Parameter value
# vec_Theta.Hat[[1]]



### MSE and Bias function (Secondary function)

"Bias.MSE" = function(Theta.hat, Theta)
{
  Bias = mean(Theta.hat - Theta)
  MSE = mean((Theta.hat - Theta)^2)
  tab = round(c(Bias, MSE), 5)
  names(tab) = c("Bias", "MSE")
  # tab = round(data.frame(Bias = Bias, MSE = MSE), 5)
  return(tab)
}


## Monte Carlo Simulation function for Xgamma
# Generates Average Bias and MSE of the estimator Theta.hat



"simple.Sim.xgamma" = function(N, n, Theta)
{
  # Generate Xgamma observations
  Xgamma_data = sapply(
    1:N,
    function(i) rxgamma(n=n, Theta = Theta)
  )
  
  # Maximum Likelihood 
  res.Theta.hat = apply(
    Xgamma_data, 
    MARGIN = 2, # function applied to columns
    function(x)
    {
      optim(par = 0.2,
            fn = loglik_Xgamma, 
            x = x, 
            method = c("Brent"), 
            lower = 0, upper = 3000,
            hessian = T
            )
    }
  )
  
  # Finding Theta.hat
  Theta.hat = sapply(
    1:N,
    function(x)res.Theta.hat[[x]]$par
  )
  
  # Bias and MSE computation
  tab = Bias.MSE(Theta.hat = Theta.hat, Theta = Theta)

  return(tab)
}


# Testing Simulation function

N = 10000
n = 20
Theta = 0.1
simple.Sim.xgamma(N=N, n=n, Theta = Theta)




# Simulation for n = 20
N = 10000
n = 20
Bias_MSE_mat = data.frame(
  t(sapply(Theta_values, 
           function(Theta)
             simple.Sim.xgamma(N=N, n=n, Theta = Theta))))



### Calling simulation function for selected values of n and Theta

n.Theta = expand.grid(Theta = Theta_values, n = n_values)
N = rep(10000, length(n.Theta))



### Simulation function for all values of n and Theta

"Sim.xgamma" = function(N, n, Theta)
{
  Bias.MSE_matrix = data.frame(
    t(
      mapply(simple.Sim.xgamma, N=N, n=n, Theta = Theta)), 
    row.names = paste("n", n, "_Theta", Theta, sep = ""))

  return(Bias.MSE_matrix)
}


## Monte Carlo Simulation function for exponential

# Simulation function for exponential distribution

"simple.Sim.exp" = function(N, n, Theta)
{
  # Generate random data from exponential distribution
  exp.Data = rexp(n = n, rate = Theta)
  
  # Maximum Likelihood estimation
  res.Theta.hat = fitdistr(
    x = exp.Data, 
    densfun = "exponential"
  )
  
  # Return MLE
  Theta.hat = res.Theta.hat$estimate
  
  # Compute Bias and MSE
  tab = Bias.MSE(Theta.hat = Theta.hat, Theta = Theta)
  
  return(tab)
}


# Testing exponential distribution simulation function

N = 10000
n = 20
Theta = 0.1
simple.Sim.exp(N=N, n=n, Theta = Theta)


### Main simulation function for exponential

"Sim.exp" = function(N, n, Theta)
{
  Bias.MSE_matrix = data.frame(
    t(
      mapply(simple.Sim.exp, N=N, n=n, Theta = Theta)), 
    row.names = paste("n", n, "_Theta", Theta, sep = ""))

  return(Bias.MSE_matrix)
}



## Simulations results - Bias and MSE 

# For Xgamma

Bias.MSE.xgamma = Sim.xgamma(N, n.Theta$n, n.Theta$Theta)

Bias.MSE.xgamma.df = data.frame(
  n = rep(n_values, each = 6),
  Theta = rep(Theta_values, times = 3),
  Bias = Bias.MSE.xgamma$Bias,
  MSE = Bias.MSE.xgamma$MSE
)

Sim.Bias.MSE.xgamma = Bias.MSE.xgamma.df %>%
  pivot_wider(names_from = n, 
              values_from = c(Bias, MSE),
              names_vary = "slowest")




# For Exponential

Bias.MSE.exp = Sim.exp(N, n.Theta$n, n.Theta$Theta)

Bias.MSE.exp.df = data.frame(
  n = rep(n_values, each = 6),
  Theta = rep(Theta_values, times = 3),
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
