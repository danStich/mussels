# Library load ----
library(R2jags)
library(tidyverse)
library(reshape)
library(lubridate)


# Data read ----
sabine <- read.csv("sabine/data/sabine.csv")
sabine_species <- read.csv("sabine/data/species_codes.csv")
sabine_species$species_2 <- gsub(" ", "", sabine_species$species_2)
individual_species <- read.csv("sabine/data/individual_species.csv")
individual_species$species_2 <- gsub(" ", "", individual_species$species_2)

# Have a look
glimpse(sabine)

# Remove individuals without a Hallprint tag
mussels <- sabine %>% 
  filter(hallprint1 != "") %>% 
  mutate(individual = hallprint1)

# Add a column for year
mussels$year <- lubridate::year(as.Date(mussels$date, format = "%m/%d/%Y"))

# Add a column for binomials
mussels <- merge(mussels, individual_species, by = "individual")

# Fix erroneous species
mussels$species_2[mussels$species_2 == "LEFR"] <- "POFR"
mussels$species_2[mussels$species_2 == "QUQU"] <- "QUAP"

# Remove individuals with missing species (there are several in data)
mussels <- mussels %>%
  filter(species_2 != "")

# Add binomials and common names
mussels <- merge(mussels, sabine_species, by = "species_2")

# Group sub-sites 
mussels$site <- "upper"
mussels$site[mussels$site_no %in% c(3, 4)] <- "lower"


# Table 2 (summary of counts) ----
sabine_species_counts <- mussels %>%
  group_by(binomial, common, site) %>% 
  summarize(n = n_distinct(individual))

sabine_species_table <- sabine_species_counts  %>%
  pivot_wider(names_from = c(site), values_from = n)

write.table(sabine_species_table, "sabine/results/species_counts.csv",
            row.names = FALSE, quote = FALSE, sep = ",")


# Capture histories ----
caps <- reshape::cast(mussels, formula = individual ~ year)

y_mat <- as.matrix(caps[ , 2:4]) 
y_mat[y_mat > 0] <- 1

# Covariates ----
covs <- mussels[!duplicated(mussels$individual), ]

# Detection Models ----
# . Shared MCMC settings ----
n_iter <- 50000
n_burnin <- 40000
chains <- 3
thin <- 5

# . Phi.P. Model ----
# .. Package the data ----
# Create vector with occasion of marking for each mussel
get.first <- function(x) min(which(x!=0))
f <- apply(y_mat, 1, get.first)

jags.data <- list(
  n_ind = nrow(y_mat),
  n_periods = ncol(y_mat),
  y = y_mat,
  f = f
)

# .. Initial values ----
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
  list(mean_phi = rnorm(1, 0, 1),
       mean_p = rnorm(1, 0, 1),
       z = known.state.cjs(y_mat))
}

# .. Fit the model ----
phi_dot_p_dot <- jags(data = jags.data, 
                      inits = inits,
                      parameters.to.save = c("mean_phi", "mean_p"),
                      model.file = "sabine/models/phi_dot_p_dot",
                      n.iter = n_iter,
                      n.burnin = n_burnin,
                      n.chains = chains,
                      n.thin = thin)

print(phi_dot_p_dot)

save(phi_dot_p_dot, file = "sabine/results/phi_dot_p_dot_50k.rda")

# . Phi.Pt Model ----
# .. Package the data ----
# Create vector with occasion of marking for each mussel
get.first <- function(x) min(which(x!=0))
f <- apply(y_mat, 1, get.first)

jags.data <- list(
  n_ind = nrow(y_mat),
  n_periods = ncol(y_mat),
  y = y_mat,
  f = f
)

# .. Initial values ----
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
  list(mean_phi = runif(1, 0, 1),
       mean_p = runif(1, 0, 1),
       z = known.state.cjs(y_mat))
}

# .. Fit the model ----
phi_dot_p_t <- jags(data = jags.data, 
                    inits = inits,
                    parameters.to.save = c("mean_phi", "mean_p", "lp"),
                    model.file = "sabine/models/phi_dot_p_t",
                    n.iter = n_iter,
                    n.burnin = n_burnin,
                    n.chains = chains,
                    n.thin = thin)

print(phi_dot_p_t)

save(phi_dot_p_t, file = "sabine/results/phi_dot_p_t_50k.rda")

# . Phi.Ps Model ----
# .. Package the data ----
# Create vector with occasion of marking for each mussel
get.first <- function(x) min(which(x!=0))
f <- apply(y_mat, 1, get.first)

jags.data <- list(
  n_ind = nrow(y_mat),
  n_periods = ncol(y_mat),
  n_species = length(unique(covs$binomial)),
  species = as.numeric(as.factor(covs$binomial)),
  y = y_mat,
  f = f
)

# .. Initial values ----
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
  list(mean_phi = runif(1, 0, 1),
       mean_p = runif(1, 0, 1),
       z = known.state.cjs(y_mat))
}

# .. Fit the model ----
phi_dot_p_s <- jags(data = jags.data, 
                    inits = inits,
                    parameters.to.save = c("mean_phi", "mean_p", "lp"),
                    model.file = "sabine/models/phi_dot_p_s",
                    n.iter = n_iter,
                    n.burnin = n_burnin,
                    n.chains = chains,
                    n.thin = thin)

print(phi_dot_p_s)

save(phi_dot_p_s, file = "sabine/results/phi_dot_p_s_50k.rda")

# . Phi.Pj Model ----
# .. Package the data ----
# Create vector with occasion of marking for each mussel
get.first <- function(x) min(which(x!=0))
f <- apply(y_mat, 1, get.first)

jags.data <- list(
  n_ind = nrow(y_mat),
  n_periods = ncol(y_mat),
  n_site = length(unique(covs$site)),
  site = as.numeric(as.factor(covs$site)),
  y = y_mat,
  f = f
)

# .. Initial values ----
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
  list(mean_phi = runif(1, 0, 1),
       mean_p = runif(1, 0, 1),
       z = known.state.cjs(y_mat))
}

# .. Fit the model ----
phi_dot_p_j <- jags(data = jags.data, 
                    inits = inits,
                    parameters.to.save = c("mean_phi", "mean_p", "lp"),
                    model.file = "sabine/models/phi_dot_p_j",
                    n.iter = n_iter,
                    n.burnin = n_burnin,
                    n.chains = chains,
                    n.thin = thin)

print(phi_dot_p_j)

save(phi_dot_p_j, file = "sabine/results/phi_dot_p_j_50k.rda")

# . Phi.Pts Model ----
# .. Package the data ----
# Create vector with occasion of marking for each mussel
get.first <- function(x) min(which(x!=0))
f <- apply(y_mat, 1, get.first)

jags.data <- list(
  n_ind = nrow(y_mat),
  n_periods = ncol(y_mat),
  n_species = length(unique(covs$binomial)),
  species = as.numeric(as.factor(covs$binomial)),
  y = y_mat,
  f = f
)

# .. Initial values ----
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
  list(mean_phi = runif(1, 0, 1),
       mean_p = runif(1, 0, 1),
       z = known.state.cjs(y_mat))
}

# .. Fit the model ----
phi_dot_p_ts <- jags(data = jags.data, 
                     inits = inits,
                     parameters.to.save = c("mean_phi", "mean_p", "lp"),
                     model.file = "sabine/models/phi_dot_p_ts",
                     n.iter = n_iter,
                     n.burnin = n_burnin,
                     n.chains = chains,
                     n.thin = thin)

print(phi_dot_p_ts)

save(phi_dot_p_ts, file = "sabine/results/phi_dot_p_ts_50k.rda")


# . Phi.Ptj Model ----
# .. Package the data ----
# Create vector with occasion of marking for each mussel
get.first <- function(x) min(which(x!=0))
f <- apply(y_mat, 1, get.first)

jags.data <- list(
  n_ind = nrow(y_mat),
  n_periods = ncol(y_mat),
  n_site = length(unique(covs$site)),
  site = as.numeric(as.factor(covs$site)),
  y = y_mat,
  f = f
)

# .. Initial values ----
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
  list(mean_phi = runif(1, 0, 1),
       mean_p = runif(1, 0, 1),
       z = known.state.cjs(y_mat))
}

# .. Fit the model ----
phi_dot_p_tj <- jags(data = jags.data, 
                     inits = inits,
                     parameters.to.save = c("mean_phi", "mean_p", "lp"),
                     model.file = "sabine/models/phi_dot_p_tj",
                     n.iter = n_iter,
                     n.burnin = n_burnin,
                     n.chains = chains,
                     n.thin = thin)

print(phi_dot_p_tj)

save(phi_dot_p_tj, file = "sabine/results/phi_dot_p_tj_50k.rda")


# . Phi.Ptj Model ----
# .. Package the data ----
# Create vector with occasion of marking for each mussel
get.first <- function(x) min(which(x!=0))
f <- apply(y_mat, 1, get.first)

jags.data <- list(
  n_ind = nrow(y_mat),
  n_periods = ncol(y_mat),
  n_site = length(unique(covs$site)),
  site = as.numeric(as.factor(covs$site)),
  n_species = length(unique(covs$binomial)),
  species = as.numeric(as.factor(covs$binomial)),  
  y = y_mat,
  f = f
)

# .. Initial values ----
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
  list(mean_phi = runif(1, 0, 1),
       mean_p = runif(1, 0, 1),
       z = known.state.cjs(y_mat))
}

# .. Fit the model ----
phi_dot_p_tsj <- jags(data = jags.data, 
                      inits = inits,
                      parameters.to.save = c("mean_phi", "mean_p", "lp"),
                      model.file = "sabine/models/phi_dot_p_tsj",
                      n.iter = n_iter,
                      n.burnin = n_burnin,
                      n.chains = chains,
                      n.thin = thin)

print(phi_dot_p_tsj)

save(phi_dot_p_tsj, file = "sabine/results/phi_dot_p_tsj_50k.rda")


# Survival models ----
# . Phi_t_p_s Model ----
# .. Package the data ----
# Create vector with occasion of marking for each mussel
get.first <- function(x) min(which(x!=0))
f <- apply(y_mat, 1, get.first)

jags.data <- list(
  n_ind = nrow(y_mat),
  n_periods = ncol(y_mat),
  n_species = length(unique(covs$binomial)),
  species = as.numeric(as.factor(covs$binomial)),
  y = y_mat,
  f = f
)

# .. Initial values ----
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
  list(mean_phi = rnorm(1, 0, 1),
       mean_p = rnorm(1, 0, 1),
       z = known.state.cjs(y_mat))
}

# .. Fit the model ----
phi_t_p_s <- jags(data = jags.data, 
                    inits = inits,
                    parameters.to.save = c("mean_phi", "mean_p", "lphi", "lp"),
                    model.file = "sabine/models/phi_t_p_s",
                    n.iter = n_iter,
                    n.burnin = n_burnin,
                    n.chains = chains,
                    n.thin = thin)

print(phi_t_p_s)

save(phi_t_p_s, file = "sabine/results/phi_t_p_s_50k.rda")

# . Phi_j_p_s Model ----
# .. Package the data ----
# Create vector with occasion of marking for each mussel
get.first <- function(x) min(which(x!=0))
f <- apply(y_mat, 1, get.first)

jags.data <- list(
  n_ind = nrow(y_mat),
  n_periods = ncol(y_mat),
  n_site = length(unique(covs$site)),
  site = as.numeric(as.factor(covs$site)),
  n_species = length(unique(covs$binomial)),
  species = as.numeric(as.factor(covs$binomial)),
  y = y_mat,
  f = f
)

# .. Initial values ----
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
  list(mean_phi = rnorm(1, 0, 1),
       mean_p = rnorm(1, 0, 1),
       z = known.state.cjs(y_mat))
}

# .. Fit the model ----
phi_j_p_s <- jags(data = jags.data, 
                    inits = inits,
                    parameters.to.save = c("mean_phi", "mean_p", "lphi"),
                    model.file = "sabine/models/phi_j_p_s",
                    n.iter = n_iter,
                    n.burnin = n_burnin,
                    n.chains = chains,
                    n.thin = thin)

print(phi_j_p_s)

save(phi_j_p_s, file = "sabine/results/phi_j_p_s_50k.rda")



# . Phi_s_p_s Model ----
# .. Package the data ----
# Create vector with occasion of marking for each mussel
get.first <- function(x) min(which(x!=0))
f <- apply(y_mat, 1, get.first)

jags.data <- list(
  n_ind = nrow(y_mat),
  n_periods = ncol(y_mat),
  n_species = length(unique(covs$binomial)),
  species = as.numeric(as.factor(covs$binomial)),
  y = y_mat,
  f = f
)

# .. Initial values ----
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
  list(mean_phi = rnorm(1, 0, 1),
       mean_p = rnorm(1, 0, 1),
       z = known.state.cjs(y_mat))
}

# .. Fit the model ----
phi_s_p_s <- jags(data = jags.data, 
                    inits = inits,
                    parameters.to.save = c("mean_phi", "mean_p", "lphi"),
                    model.file = "sabine/models/phi_s_p_s",
                    n.iter = n_iter,
                    n.burnin = n_burnin,
                    n.chains = chains,
                    n.thin = thin)

print(phi_s_p_s)

save(phi_s_p_s, file = "sabine/results/phi_s_p_s_50k.rda")


# . Phi_tj_p_s Model ----
# .. Package the data ----
# Create vector with occasion of marking for each mussel
get.first <- function(x) min(which(x!=0))
f <- apply(y_mat, 1, get.first)

jags.data <- list(
  n_ind = nrow(y_mat),
  n_periods = ncol(y_mat),
  n_site = length(unique(covs$site)),
  site = as.numeric(as.factor(covs$site)),
  n_species = length(unique(covs$binomial)),
  species = as.numeric(as.factor(covs$binomial)),
  y = y_mat,
  f = f
)

# .. Initial values ----
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
  list(mean_phi = rnorm(1, 0, 1),
       mean_p = rnorm(1, 0, 1),
       z = known.state.cjs(y_mat))
}

# .. Fit the model ----
phi_tj_p_s <- jags(data = jags.data, 
                     inits = inits,
                     parameters.to.save = c("mean_phi", "mean_p", "lphi", "lp"),
                     model.file = "sabine/models/phi_tj_p_s",
                     n.iter = n_iter,
                     n.burnin = n_burnin,
                     n.chains = chains,
                     n.thin = thin)

print(phi_tj_p_s)

save(phi_tj_p_s, file = "sabine/results/phi_tj_p_s_50k.rda")


# . Phi_ts_p_s Model ----
# .. Package the data ----
# Create vector with occasion of marking for each mussel
get.first <- function(x) min(which(x!=0))
f <- apply(y_mat, 1, get.first)

jags.data <- list(
  n_ind = nrow(y_mat),
  n_periods = ncol(y_mat),
  n_species = length(unique(covs$binomial)),
  species = as.numeric(as.factor(covs$binomial)),
  y = y_mat,
  f = f
)

# .. Initial values ----
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
  list(mean_phi = rnorm(1, 0, 1),
       mean_p = rnorm(1, 0, 1),
       z = known.state.cjs(y_mat))
}

# .. Fit the model ----
phi_ts_p_s <- jags(data = jags.data, 
                     inits = inits,
                     parameters.to.save = c("mean_phi", "mean_p", "lphi", "lp"),
                     model.file = "sabine/models/phi_ts_p_s",
                     n.iter = n_iter,
                     n.burnin = n_burnin,
                     n.chains = chains,
                     n.thin = thin)

print(phi_ts_p_s)

save(phi_ts_p_s, file = "sabine/results/phi_ts_p_s_50k.rda")


# . Phi_js_p_s Model ----
# .. Package the data ----
# Create vector with occasion of marking for each mussel
get.first <- function(x) min(which(x!=0))
f <- apply(y_mat, 1, get.first)

jags.data <- list(
  n_ind = nrow(y_mat),
  n_periods = ncol(y_mat),
  n_species = length(unique(covs$binomial)),
  species = as.numeric(as.factor(covs$binomial)),
  n_site = length(unique(covs$site)),
  site = as.numeric(as.factor(covs$site)),  
  y = y_mat,
  f = f
)

# .. Initial values ----
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
  list(mean_phi = rnorm(1, 0, 1),
       mean_p = rnorm(1, 0, 1),
       z = known.state.cjs(y_mat))
}

# .. Fit the model ----
phi_js_p_s <- jags(data = jags.data, 
                     inits = inits,
                     parameters.to.save = c("mean_phi", "mean_p", "lphi", "lp"),
                     model.file = "sabine/models/phi_js_p_s",
                     n.iter = n_iter,
                     n.burnin = n_burnin,
                     n.chains = chains,
                     n.thin = thin)

print(phi_js_p_s)

save(phi_js_p_s, file = "sabine/results/phi_js_p_s_50k.rda")


# . Phi_tjs_p_s Model ----
# .. Package the data ----
# Create vector with occasion of marking for each mussel
get.first <- function(x) min(which(x!=0))
f <- apply(y_mat, 1, get.first)

jags.data <- list(
  n_ind = nrow(y_mat),
  n_periods = ncol(y_mat),
  n_species = length(unique(covs$binomial)),
  species = as.numeric(as.factor(covs$binomial)),
  n_site = length(unique(covs$site)),
  site = as.numeric(as.factor(covs$site)),  
  y = y_mat,
  f = f
)

# .. Initial values ----
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
  list(mean_phi = rnorm(1, 0, 1),
       mean_p = rnorm(1, 0, 1),
       z = known.state.cjs(y_mat))
}

# .. Fit the model ----
phi_tjs_p_s <- jags(data = jags.data, 
                      inits = inits,
                      parameters.to.save = c("mean_phi", "mean_p", "lphi", "lp"),
                      model.file = "sabine/models/phi_tjs_p_s",
                      n.iter = n_iter,
                      n.burnin = n_burnin,
                      n.chains = chains,
                      n.thin = thin)

print(phi_tjs_p_s)

save(phi_tjs_p_s, file = "sabine/results/phi_tjs_p_s_50k.rda")


