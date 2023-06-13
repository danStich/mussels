# Library load ----
library(tidyverse)

# Load results ----
load(file = "simulations/plri_15k_iter.rda")

# Summarize results ----
plri_plotter <- res %>% 
  group_by(n_ind, n_periods) %>% 
  filter(n_periods <= 30) %>%
  summarize(
    error_mean = mean(obs_phi - phi),
    error_lwr = quantile(obs_phi - phi, 0.025),
    error_upr = quantile(obs_phi - phi, 0.975),
    precision_mean = mean(obs_phi_sd),
    precision_lwr = quantile(obs_phi_sd, 0.025),
    precision_upr = quantile(obs_phi_sd, 0.975),  
    mse = (sum(obs_phi - phi)^2)/n()) %>% 
  mutate(Species = "Pleurobema riddellii")


# Load results ----
load(file = "simulations/poam_15k_iter.rda")

# Summarize results ----
poam_plotter <- res %>% 
  group_by(n_ind, n_periods) %>% 
  filter(n_periods <= 30) %>%
  summarize(
    error_mean = mean(obs_phi - phi),
    error_lwr = quantile(obs_phi - phi, 0.025),
    error_upr = quantile(obs_phi - phi, 0.975),
    precision_mean = mean(obs_phi_sd),
    precision_lwr = quantile(obs_phi_sd, 0.025),
    precision_upr = quantile(obs_phi_sd, 0.975),  
    mse = (sum(obs_phi - phi)^2)/n()) %>% 
  mutate(Species = "Potamilus amphichaenus")

plotter <- rbind(plri_plotter, poam_plotter)


# Line graph for bias
ggplot(plotter, aes(x = n_ind, y = error_mean,
                    color = factor(n_periods), fill = factor(n_periods))) +
  geom_hline(yintercept = 0, lty=2) +
  scale_y_continuous(limits = c(-.2, .2)) +
  geom_line() +
  geom_ribbon(
    aes(xmax = n_ind, ymin = error_lwr, ymax = error_upr,
        color = NULL),
    alpha = 0.15) +
  ylab(expression(paste("Error (", hat(phi) - phi, ")"))) +
  xlab("Sample size") +
  labs(color = "N periods", fill = "N periods") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  facet_wrap(~Species)


# Line graph for precision
ggplot(plotter, aes(x = n_ind, y = precision_mean,
                    color = factor(n_periods), fill = factor(n_periods))) +
  scale_y_continuous(limits = c(0, .1)) +
  geom_line() +
  geom_ribbon(
    aes(xmax = n_ind, ymin = precision_lwr, ymax = precision_upr,
        color = NULL),
    alpha = 0.15) +
  ylab(expression(paste("Standard deviation of "," " ,hat(phi)))) +
  xlab("Sample size") +
  labs(color = "N periods", fill = "N periods") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  facet_wrap(~Species)

# Line graph for mse
ggplot(plotter, aes(x = n_ind, y = sqrt(mse),
                    color = factor(n_periods), fill = factor(n_periods))) +
  geom_line() +
  ylab("Mean squared error") +
  xlab("Sample size") +
  labs(color = "N periods", fill = "N periods") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  facet_wrap(~Species)


