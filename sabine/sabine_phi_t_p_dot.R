# Library load ----
library(R2jags)
library(tidyverse)
library(reshape)
library(lubridate)

# Data read ----
sabine <- read.csv("data/sabine.csv", stringsAsFactors = FALSE)

# Have a look
glimpse(sabine)

# Remove individuals without a Hallprint tag
mussels <- sabine %>% 
  filter(hallprint1 != "") %>% 
  mutate(individual = paste0(hallprint1, hallprint2),
         date = as.Date(date, format = "%m/%d/%Y"),
         year = year(date)
  )

# Get species for each individual
ind_spp <- data.frame(
  unique(cbind(mussels$individual[mussels$species != ""], 
               mussels$species[mussels$species != ""])))
names(ind_spp) <- c("individual", "species_2")

# Merge back with raw data to get species in "mussels"
mussels <- merge(mussels, ind_spp, by = "individual")


# Combine dates into samples based on original attempt at 
# individual capture histories
# mussels$sample <- "sample_1"
# mussels$sample[mussels$date %in% c("8/2/2022")] <- "sample_2"
# mussels$sample[mussels$date %in% c("8/3/2022", "8/4/2022")] <- "sample_3"

# melted <- melt(mussels, id = c("individual", "site", "year"))

caps <- reshape::cast(mussels, formula = individual ~ year)


y_mat <- as.matrix(caps[ , 2:4])  
y_mat[y_mat > 0] <- 1

# Phi_t_P. Model ----
# . Package the data ----
# Create vector with occasion of marking for each mussel
get.first <- function(x) min(which(x!=0))
f <- apply(y_mat, 1, get.first)

jags.data <- list(
  n_ind = nrow(y_mat),
  n_periods = ncol(y_mat),
  y = y_mat,
  f = f
)

# . MCMC settings ----
n_iter <- 25000
n_burnin <- 1500
chains <- 3
thin <- 3

# . Initial values ----
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

# . Fit the model ----
phi_t_p_dot <- jags(data = jags.data, 
                      inits = inits,
                      parameters.to.save = c("mean_phi", "mean_p","phi", "p"),
                      model.file = "models/phi_t_p_dot"
)

print(phi_t_p_dot)

save(phi_t_p_dot, file = "results/phi_dot_p_dot_1.rda")
