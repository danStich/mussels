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

# Model ----
# . Phi_js_p_tjs----
# .. Package the data ----
# Drop any missing rows from ch and covariate matrix
# covs <- covs[rowSums(y_mat) > 0 & y_mat[ ,1] == 1 ,]
# y_mat <- y_mat[rowSums(y_mat) > 0 & y_mat[ ,1] == 1 ,]


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
  f = f)

# .. MCMC settings ----
n_iter <- 5000
n_burnin <- 4000
chains <- 3
thin <- 5

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
# phi_ts_p_s <- jags(data = jags.data, 
#                      inits = inits,
#                      parameters.to.save = c("mean_phi", "mean_p", "lphi", "lp"),
#                      model.file = "sabine/models/phi_ts_p_s",
#                      n.iter = n_iter,
#                      n.burnin = n_burnin,
#                      n.chains = chains,
#                      n.thin = thin)
# 
# print(phi_ts_p_s)
# 

# Results----
load("sabine/results/phi_ts_p_s_50k.rda")

# phi_ts_p_s ----
# . Survival estimates ----
names(phi_ts_p_s$BUGSoutput$sims.list)

# Get logit-scale survival estimates
lphi <-  phi_ts_p_s$BUGSoutput$sims.list$lphi
lphi <- melt(lphi)
names(lphi) <- c("iteration", "time", "species", "estimate")
phi <- filter(lphi, time ==1)
phi$estimate <- boot::inv.logit(phi$estimate)
phi$fit <-phi$estimate

# Change species to factor
phi$common <- levels(factor(covs$common))[phi$species]
phi$binomial <- levels(factor(covs$binomial))[phi$species]

# Calculate descriptive stats for each species
phi_ests <- phi %>% 
  group_by(binomial, common) %>% 
  summarize(
    fit = median(estimate),
    sd = sd(estimate),
    lwr = quantile(estimate, 0.025),
    upr = quantile(estimate, 0.975),
    q1 = quantile(estimate, 0.25), 
    q2 = quantile(estimate, 0.75)) %>% 
  mutate(parameter = "phi")

# Write the descriptive statistics to a table
write.table(phi_ests, "sabine/results/phi_ests.csv",
            row.names = FALSE, quote = FALSE, sep = ",")

# Calculate the overall mean and descriptive statistics
phi_mean <- phi %>% summarize(fit = mean(estimate))

# Here they are: 
mean(phi$fit)
quantile(phi$fit, c(0.025, 0.975))

# . Survival plot ----
ggplot(phi, aes(x = fit, y = common)) +
  ggridges:::geom_density_ridges(scale = 1, size = 0.3, alpha = 0.20, 
                                 color = NA) +
  geom_point(data = phi_ests, position = position_nudge(y = 0.1), size = 2) +
  geom_linerange(data = phi_ests, 
                aes(xmin = lwr, xmax = upr),
                position = position_nudge(y = 0.1),
                linewidth = .5) +
  geom_linerange(data = phi_ests, 
                 aes(xmin = q1, xmax = q2),
                 position = position_nudge(y = 0.1),
                 linewidth = 1) +
  ylab("Common name") +
  xlab(expression(paste("Apparent survival (", phi, ")"))) + 
  scale_y_discrete(limits = rev, expand = expansion(mult = c(0.01, .06))) +
  geom_vline(data = phi_mean, aes(xintercept = fit), linetype = 2) +
  theme_light()


# . Detection estimates ----
lp <- phi_ts_p_s$BUGSoutput$sims.list$lp  
p <- melt(lp)
p$value <- boot::inv.logit(p$value)
names(p) <- c("iteration", "species", "estimate")
p$fit <- p$estimate

# Change species to factor
p$binomial <- levels(factor(covs$binomial))[p$species]
p$common <- levels(factor(covs$common))[p$species]

# Calculate descriptive stats for each species
p_ests <- p %>% 
  group_by(binomial, common) %>% 
  summarize(
    fit = median(estimate),
    sd = sd(estimate),
    lwr = quantile(estimate, 0.025),
    upr = quantile(estimate, 0.975),
    q1 = quantile(estimate, 0.25), 
    q2 = quantile(estimate, 0.75)) %>% 
  mutate(parameter = "p")

# Write the descriptive statistics to a table
write.table(p_ests, "sabine/results/p_ests.csv",
            row.names = FALSE, quote = FALSE, sep = ",")

# Calculate the overall mean and some descriptive stats
p_mean <- p %>% summarize(fit = mean(estimate))

# Here they are: 
mean(p$fit)
quantile(p$fit, c(0.025, 0.975))

# . Detection plot ----
ggplot(p, aes(x = fit, y = common)) +
  ggridges:::geom_density_ridges(scale = 1, size = 0.3, alpha = 0.20, 
                                 color = NA) +
  geom_point(data = p_ests, position = position_nudge(y = 0.1), size = 2) +
  geom_linerange(data = p_ests, 
                 aes(xmin = lwr, xmax = upr),
                 position = position_nudge(y = 0.1),
                 linewidth = .5) +
  geom_linerange(data = p_ests, 
                 aes(xmin = q1, xmax = q2),
                 position = position_nudge(y = 0.1),
                 linewidth = 1) +
  ylab("Common name") +
  xlab("Detection probability (p)") + 
  scale_y_discrete(limits = rev, expand = expansion(mult = c(0.01, .06))) +
  geom_vline(data = p_mean, aes(xintercept = fit), linetype = 2) +
  theme_light()

