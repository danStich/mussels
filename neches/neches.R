# Library load ----
library(R2jags)
library(tidyverse)
library(reshape)
library(lubridate)


# Data read ----
neches <- read.csv("neches/data/neches_20230522.csv")
neches_species <- read.csv("neches/data/species_codes.csv")
individual_species <- read.csv("neches/data/individual_species.csv")

# Have a look
glimpse(neches)

# Make a copy
mussels <- neches

# Separate initial marking event at both sites
mussels$event[mussels$date %in% c("10/5/2021", "10/6/2021")] <- 0

# Tag cleanup
mussels$hallprint1 <- gsub(" ", "", mussels$hallprint1)

# Remove individuals without a Hallprint tag
mussels <- mussels %>% 
  filter(!hallprint1 %in% c('', ' ', '-')) %>% 
  mutate(individual = hallprint1) %>% 
  arrange(individual, event, species)

# Merge back with raw data to get species in "mussels"
mussels <- merge(mussels, individual_species, by = "individual")

# Fix erroneous species
mussels$species_2[mussels$species_2 == "CYTA"] <- "GLRO"
mussels$species_2[mussels$species_2 == "ELSL"] <- "PLDO"
mussels$species_2[mussels$species_2 == "LEFR"] <- "POFR"

# Remove individuals with missing species (there are several in data)
mussels <- mussels %>% 
  filter(species_2 != "")

# Add a column for binomials - note that "TOTE" is dropped because it was
# never tagged or recaptured (there was only one observation)
mussels <- merge(mussels, neches_species)
mussels$year <- lubridate::year(as.Date(mussels$date, format = "%m/%d/%Y"))

# Table 1 (summary of counts) ----
neches_species_counts <- mussels %>%
  group_by(binomial, common, site) %>% 
  summarize(n = n_distinct(individual))

neches_species_table <- neches_species_counts %>%
  pivot_wider(names_from = c(site), values_from = n)

write.table(neches_species_table, "neches/results/species_counts.csv",
            row.names = FALSE, quote = FALSE, sep = ",")


# Capture histories ----
caps <- reshape::cast(mussels, formula = individual ~ event)

y_mat <- as.matrix(caps[ , 2:5]) 
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
                model.file = "neches/models/phi_dot_p_dot",
                n.iter = n_iter,
                n.burnin = n_burnin,
                n.chains = chains,
                n.thin = thin)

print(phi_dot_p_dot)

save(phi_dot_p_dot, file = "neches/results/phi_dot_p_dot_50k.rda")

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
                    model.file = "neches/models/phi_dot_p_t",
                    n.iter = n_iter,
                    n.burnin = n_burnin,
                    n.chains = chains,
                    n.thin = thin)

print(phi_dot_p_t)

save(phi_dot_p_t, file = "neches/results/phi_dot_p_t_50k.rda")

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
                    model.file = "neches/models/phi_dot_p_s",
                    n.iter = n_iter,
                    n.burnin = n_burnin,
                    n.chains = chains,
                    n.thin = thin)

print(phi_dot_p_s)

save(phi_dot_p_s, file = "neches/results/phi_dot_p_s_50k.rda")

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
                    model.file = "neches/models/phi_dot_p_j",
                    n.iter = n_iter,
                    n.burnin = n_burnin,
                    n.chains = chains,
                    n.thin = thin)

print(phi_dot_p_j)

save(phi_dot_p_j, file = "neches/results/phi_dot_p_j_50k.rda")

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
                    model.file = "neches/models/phi_dot_p_ts",
                    n.iter = n_iter,
                    n.burnin = n_burnin,
                    n.chains = chains,
                    n.thin = thin)

print(phi_dot_p_ts)

save(phi_dot_p_ts, file = "neches/results/phi_dot_p_ts_50k.rda")


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
                     model.file = "neches/models/phi_dot_p_tj",
                     n.iter = n_iter,
                     n.burnin = n_burnin,
                     n.chains = chains,
                     n.thin = thin)

print(phi_dot_p_tj)

save(phi_dot_p_tj, file = "neches/results/phi_dot_p_tj_50k.rda")


# . Phi.Ptsj Model ----
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
                     model.file = "neches/models/phi_dot_p_tsj",
                     n.iter = n_iter,
                     n.burnin = n_burnin,
                     n.chains = chains,
                     n.thin = thin)

print(phi_dot_p_tsj)

save(phi_dot_p_tsj, file = "neches/results/phi_dot_p_tsj_50k.rda")


# . Detection model selection ----
det_mods <- data.frame(
  Model = c("Phi.p.", "Phi.pt", "Phi.ps", "Phi.pts", "Phi.pj", "Phi.ptj",
            "Phi.ptsj"),
  DIC = c(phi_dot_p_dot$BUGSoutput$DIC,
          phi_dot_p_t$BUGSoutput$DIC,
          phi_dot_p_s$BUGSoutput$DIC,
          phi_dot_p_ts$BUGSoutput$DIC,
          phi_dot_p_j$BUGSoutput$DIC,
          phi_dot_p_tj$BUGSoutput$DIC,
          phi_dot_p_tsj$BUGSoutput$DIC))

det_mods[order(det_mods$DIC), ]

# N    Model      DIC
# 1   Phi.p. 6024.987
# 5   Phi.pj 6104.997
# 3   Phi.ps 7606.663
# 7 Phi.ptsj 8009.046
# 2   Phi.pt 8514.412
# 4  Phi.pts 8547.261
# 6  Phi.ptj 8605.194


# Survival models ----
# . Phi_t_p_dot Model ----
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
phi_t_p_dot <- jags(data = jags.data, 
                   inits = inits,
                   parameters.to.save = c("mean_phi", "mean_p", "lphi", "lp"),
                   model.file = "neches/models/phi_t_p_dot",
                   n.iter = n_iter,
                   n.burnin = n_burnin,
                   n.chains = chains,
                   n.thin = thin)

print(phi_t_p_dot)

save(phi_t_p_dot, file = "neches/results/phi_t_p_dot_50k.rda")

# . Phi_j_p_dot Model ----
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
  list(mean_phi = rnorm(1, 0, 1),
       mean_p = rnorm(1, 0, 1),
       z = known.state.cjs(y_mat))
}

# .. Fit the model ----
phi_j_p_dot <- jags(data = jags.data, 
                    inits = inits,
                    parameters.to.save = c("mean_phi", "mean_p", "lphi"),
                    model.file = "neches/models/phi_j_p_dot",
                    n.iter = n_iter,
                    n.burnin = n_burnin,
                    n.chains = chains,
                    n.thin = thin)

print(phi_j_p_dot)

save(phi_j_p_dot, file = "neches/results/phi_j_p_dot_50k.rda")



# . Phi_s_p_dot Model ----
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
phi_s_p_dot <- jags(data = jags.data, 
                    inits = inits,
                    parameters.to.save = c("mean_phi", "mean_p", "lphi"),
                    model.file = "neches/models/phi_s_p_dot",
                    n.iter = n_iter,
                    n.burnin = n_burnin,
                    n.chains = chains,
                    n.thin = thin)

print(phi_s_p_dot)

save(phi_s_p_dot, file = "neches/results/phi_s_p_dot_50k.rda")


# . Phi_tj_p_dot Model ----
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
  list(mean_phi = rnorm(1, 0, 1),
       mean_p = rnorm(1, 0, 1),
       z = known.state.cjs(y_mat))
}

# .. Fit the model ----
phi_tj_p_dot <- jags(data = jags.data, 
                    inits = inits,
                    parameters.to.save = c("mean_phi", "mean_p", "lphi", "lp"),
                    model.file = "neches/models/phi_tj_p_dot",
                    n.iter = n_iter,
                    n.burnin = n_burnin,
                    n.chains = chains,
                    n.thin = thin)

print(phi_tj_p_dot)

save(phi_tj_p_dot, file = "neches/results/phi_tj_p_dot_50k.rda")


# . Phi_ts_p_dot Model ----
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
phi_ts_p_dot <- jags(data = jags.data, 
                     inits = inits,
                     parameters.to.save = c("mean_phi", "mean_p", "lphi", "lp"),
                     model.file = "neches/models/phi_ts_p_dot",
                     n.iter = n_iter,
                     n.burnin = n_burnin,
                     n.chains = chains,
                     n.thin = thin)

print(phi_ts_p_dot)

save(phi_ts_p_dot, file = "neches/results/phi_ts_p_dot_50k.rda")


# . Phi_js_p_dot Model ----
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
phi_js_p_dot <- jags(data = jags.data, 
                     inits = inits,
                     parameters.to.save = c("mean_phi", "mean_p", "lphi", "lp"),
                     model.file = "neches/models/phi_js_p_dot",
                     n.iter = n_iter,
                     n.burnin = n_burnin,
                     n.chains = chains,
                     n.thin = thin)

print(phi_js_p_dot)

save(phi_js_p_dot, file = "neches/results/phi_js_p_dot_50k.rda")


# . Phi_tjs_p_dot Model ----
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
phi_tjs_p_dot <- jags(data = jags.data, 
                     inits = inits,
                     parameters.to.save = c("mean_phi", "mean_p", "lphi", "lp"),
                     model.file = "neches/models/phi_tjs_p_dot",
                     n.iter = n_iter,
                     n.burnin = n_burnin,
                     n.chains = chains,
                     n.thin = thin)

print(phi_tjs_p_dot)

save(phi_tjs_p_dot, file = "neches/results/phi_tjs_p_dot_50k.rda")

