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
# phi_tsj_p_dot <- jags(data = jags.data, 
#                      inits = inits,
#                      parameters.to.save = c("mean_phi", "mean_p", "lphi", "lp"),
#                      model.file = "neches/models/phi_tsj_p_dot",
#                      n.iter = n_iter,
#                      n.burnin = n_burnin,
#                      n.chains = chains,
#                      n.thin = thin)
# 
# print(phi_tsj_p_dot)


# Results----
# . Load result ----
load("neches/results/phi_tsj_p_dot_50k.rda")

# . Phi_tsj_p_dot ----
# .. Survival ----
names(phi_tsj_p_dot$BUGSoutput$sims.list)

# ... Estimates ----
# Get logit-scale survival estimates
lphi <-  phi_tsj_p_dot$BUGSoutput$sims.list$lphi
lphi <- melt(lphi)
names(lphi) <- c("iteration", "time", "site", "species", "estimate")
phi <- filter(lphi, time ==1)
phi$estimate <- boot::inv.logit(phi$estimate)
phi$fit <- phi$estimate


# Change species to factor
phi$common <- levels(factor(covs$common))[phi$species]
phi$binomial <- levels(factor(covs$binomial))[phi$species]

phi$site <- levels(factor(covs$site))[phi$site]


# Calculate descriptive stats for each species
phi_ests <- phi %>% 
  group_by(binomial, common, site) %>% 
  summarize(
    fit = median(estimate),
    sd = sd(estimate),
    lwr = quantile(estimate, 0.025),
    upr = quantile(estimate, 0.975),
    q1 = quantile(estimate, 0.25), 
    q2 = quantile(estimate, 0.75)) %>% 
  mutate(parameter = "phi")

# Write the descriptive statistics to a table
write.table(phi_ests, "neches/results/phi_ests.csv",
            row.names = FALSE, quote = FALSE, sep = ",")

# Calculate the overall mean and descriptive statistics
phi_mean <- phi %>% 
  group_by(site) %>% 
  summarize(fit = median(estimate),
            lwr = quantile(estimate, 0.025),
            upr = quantile(estimate, 0.975))

# Here is the overall:
median(phi$fit)
quantile(phi$fit, c(0.025, 0.975))

# ... Plot ----
# Create an index to determine whether to plot the 
# point estimates and quantiles for species that were not detected
# at a given site
# Do the same thing for the full posteriors
include_df <- neches_species_table %>% 
  pivot_longer(cols = c(`Lower Neches`, Rockland), names_to = "site") %>% 
  mutate(value = ifelse(is.na(value), 0, 1)) %>% 
  data.frame()

# Drop rows from estimates with zero mussels in a site, since these
# are just the means from the other site
phi_ests <- phi_ests[phi_ests$fit * include_df$value > 0, ]

phi2 <- phi
phi3 <- merge(phi2, include_df, by = c("binomial", "site")) %>% 
  filter(value > 0) %>% 
  select(-common.y) %>% 
  dplyr::rename(common = common.x)


# Plot the full posteriors as a density ridge overlayed by
# point estimates and quantiles
ggplot(phi3, aes(x = fit, y = common)) +
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
  theme_light() +
  facet_wrap(~site)

# .. Detection probability ----
# ... Estimates ----
# Get logit-scale survival estimates
lp <-  phi_tsj_p_dot$BUGSoutput$sims.list$lp
p <- melt(lp)
names(p) <- c("iteration", "time", "site", "species", "estimate")
p$estimate <- boot::inv.logit(p$estimate)
p$fit <- p$estimate

# Change species to factor
p$common <- levels(factor(covs$common))[p$species]
p$binomial <- levels(factor(covs$binomial))[p$species]

p$site <- levels(factor(covs$site))[p$site]

# Calculate the overall mean and descriptive statistics
p_ests <- p %>% 
  summarize(fit = median(estimate),
            lwr = quantile(estimate, 0.025),
            upr = quantile(estimate, 0.975),
            q1 = quantile(estimate, 0.25),
            q2 = quantile(estimate, 0.75))

# ... Plot ----
ggplot(p, aes(x = fit, y = 1)) +
  ggridges:::geom_density_ridges(scale = 1, size = 0.3, alpha = 0.20, 
                                 color = NA) +
  geom_point(data = p_ests, position = position_nudge(y = 1.5), size = 4) +
  geom_linerange(data = p_ests, 
                 aes(xmin = lwr, xmax = upr),
                 position = position_nudge(y = 1.5),
                 linewidth = 1) +
  geom_linerange(data = p_ests, 
                 aes(xmin = q1, xmax = q2),
                 position = position_nudge(y = 1.5),
                 linewidth = 2) +
  scale_y_continuous(expand = c(0, 1)) +
  ylab("Density") +
  xlab("Detection probability (p)") + 
  theme_light()
