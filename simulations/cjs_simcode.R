# Library load ----
library(R2jags)

library(snowfall)

# Set up a parallel back-end ----
# Initialize snowfall
sfInit(parallel=TRUE, cpus=7, type="SOCK")

# Wrapper function ----
wrapper <- function(x){
  
  # . Define known values ----
  # How to simulate 1 and 0 data
  # from survival and detection probabilities
  
  # Simulate 1 or 0 for each of n individuals in
  # 10 time periods (observation events):
  
  # We need to define some variables:
  # Number of observations
  periods <- seq(30, 50, 5) #seq(5, 50, 5)
  n_periods <- sample(periods, 1, replace = TRUE)
  
  # Number of individuals
  inds <- c(seq(10, 100, 10), 200, 300, 400, 500)
  n_ind <- sample(inds, 1, replace = TRUE)
  
  # Known survival probability
  # A single survival rate for all individuals and periods
  phi <- 0.83
  
  # Known detection probability
  # Also same across individuals and time
  p <- 0.65
  
  # . Create random samples based on variables above ----
  # .. State matrix ----
  # Create blank matrix to hold survival outcome ("state")
  z_mat <- matrix(data = NA, nrow = n_ind, ncol = n_periods)
  
  # All individuals released at time 1 in this example,
  # so we add a "1" to the first column in all rows
  z_mat[, 1] <- 1
  
  # For each row, and each time period from 2 until the
  # end, state is the outcome (1 or 0) of a random binomial
  # draw with size = state in previous time and probability of
  # success = survival (phi)
  for(i in 1:nrow(z_mat)){
    for(t in 2:ncol(z_mat)){
      z_mat[i, t] <- rbinom(n = 1, size = z_mat[i, t - 1], prob = phi)
    } 
  }
  
  # .. Observation matrix ----
  y_mat <- z_mat
  for(i in 1:nrow(y_mat)){
    for(t in 2:ncol(y_mat)){
      y_mat[i, t] <- rbinom(n = 1, size = y_mat[i, t], prob = p)
    } 
  }
  
  
  # . JAGS Model ----
  # . Cormack-Jolly-Seber ----
  model_string = "
  model{
  
    # Likelihood
    for(i in 1:n_ind){
    
      # Define latent state at first capture. For this example
      # they are all first captured in time (column) 1. This was 
      # missing in our first attempt. We will re-write to use the
      # first occasion once we are up and running
      z[i, 1] <- 1

      for(t in 2:n_periods){
         # Z (true state) as outcome of observed state in
         # previous time period and survival
         z[i, t] ~ dbern(state[i, t])
         state[i, t] <- phi * z[i, t - 1]
         
         # y (observed state - data) conditional on true state
         y[i, t] ~ dbern(observed[i, t])
         observed[i, t] <- p * z[i, t]
      }
    }
    
    phi ~ dbeta(1, 1)
    p ~ dbeta(1, 1)
  
  }"
  
  # . Run the model ----
  # .. Package the data ----
  jags.data <- list(
    n_ind = n_ind,
    n_periods = n_periods,
    y = y_mat
  )
  
  # .. MCMC settings ----
  n_iter <- 100000
  n_burnin <- 80000
  chains <- 3
  thin <- 3
  
  # .. Initial values ----
  # This was the second problem
  # In JAGS we have to give good initial values for the latent state z. 
  # At all occasions when an individual was observed, its state is z = 1 
  # for sure. In addition, if an individual was not observed at an 
  # occasion, but was alive for sure, because it was observed before 
  # and thereafter (i.e. has a capture history of e.g. {101} or {10001}), 
  # then we know that the individual was alive at all of these occasions, 
  # and thus z = 1. Therefore, we should provide initial values of z = 1 
  # at these positions as well. The following function provides such 
  # initial values from the observed capture histories:
  known.state.cjs <- function(ch){
    state <- ch
    for (i in 1:dim(ch)[1]){
      n1 <- min(which(ch[i,]==1))
      n2 <- max(which(ch[i,]==1))
      state[i,n1:n2] <- 1
      state[i,n1] <- NA
    }
    state[state==0] <- NA
    return(state)
  }
  
  
  # This line creates a function that creates "initial values"
  # for the jags model. For phi and p we are drawing random 
  # values from a uniform on (0, 1). For z (true state) we 
  # use the gnarly function defined above from Kery and Schaub
  inits <- function(){
    list(phi = runif(1, 0, 1),
         p = runif(1, 0, 1),
         z = known.state.cjs(y_mat))
  }
  
  # .. Fit the model ----
  cjs_out <- jags(data = jags.data, 
                  inits = inits,
                  parameters.to.save = c("phi", "p"),
                  model.file = textConnection(model_string)
  )
  
  # . Results ----
  # All the good stuff (posteriors) is in the 
  # cjs_out object
  
  # .. Extract the posterior estimates ----
  # Descriptive statistics for survival
  obs_phi <- mean(cjs_out$BUGSoutput$sims.list$phi)
  obs_phi_sd <- sd(cjs_out$BUGSoutput$sims.list$phi)
  
  # Descriptive statistics for detection
  obs_p <- mean(cjs_out$BUGSoutput$sims.list$p)
  obs_p_sd <- sd(cjs_out$BUGSoutput$sims.list$p)
  
  # .. Write true and estimated values to list ----
  sim <- list(
    phi = phi, 
    obs_phi = obs_phi, 
    obs_phi_sd = obs_phi_sd, 
    p = p, 
    obs_p = obs_p, 
    obs_p_sd = obs_p_sd,
    n_periods = n_periods,
    n_ind = n_ind
  )
  
}

# Simulation settings ----
# Necessary libraries
sfLibrary(R2jags)

# Number of simulations
n_iterations <- 15000

# Start timer
start_time <- Sys.time()

# Run simulation ----
# Run the simulation n_iterations times
result <- sfLapply(1:n_iterations, wrapper) 

# Stop snowfall
sfStop()

# Stop timer
total_time <- Sys.time()- start_time
total_time

# Process results ----
# Extract results list from output list
out <- lapply(result, data.frame)
res <- do.call(rbind, out)

# Write the results to a file ----
save(res, file = "poam_100k_iter_.rda")
