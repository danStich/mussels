# Library load ----
library(tidyverse)

# load in data ----
load( file = "simulations_results/plri_15k_iter.rda")

plotter <- res %>% 
  group_by(n_ind, n_periods) %>% 
  # filter(n_periods == 5) %>% 
  summarize(
    error_mean = mean(obs_phi - phi),
    error_lwr = quantile(obs_phi - phi, 0.025),
    error_upr = quantile(obs_phi - phi, 0.975),
    precision_mean = mean(obs_phi_sd),
    precision_lwr = quantile(obs_phi_sd, 0.025),
    precision_upr = quantile(obs_phi_sd, 0.975),  
    mse = (sum(obs_phi - phi)^2)/n()
  )

jpeg("plri_accuracy_fig.jpeg",
     width = 2000,
     height = 1500,
     res =  300)

# Line graph for bias
ggplot(plotter, aes(x = n_ind, y = error_mean,
                    color = factor(n_periods), fill = factor(n_periods)
)) +
  geom_hline(yintercept = 0, lty=2) +
  # scale_y_continuous(limits = c(-.1, .1)) +
  geom_line() +
  geom_ribbon(
    aes(xmax = n_ind, ymin = error_lwr, ymax = error_upr,
        color = NULL),
    alpha = 0.15) +
  ylab(expression(paste("Error (", hat(phi) - phi, ")"))) +
  xlab("Sample size") +
  labs(color = "N periods", fill = "N periods") +
  theme(axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 14))

dev.off()

jpeg("plri_precision_fig.jpeg",
     width = 2000,
     height = 1500,
     res =  300)


# Line graph for precision
ggplot(plotter, aes(x = n_ind, y = precision_mean,
                    color = factor(n_periods), fill = factor(n_periods)
)) +
  # geom_hline(yintercept = 0, lty=2) +
  scale_y_continuous(limits = c(0, .1)) +
  geom_line() +
  geom_ribbon(
    aes(xmax = n_ind, ymin = precision_lwr, ymax = precision_upr,
        color = NULL),
    alpha = 0.15) +
  ylab(expression(paste("Standard deviation of "," " ,hat(phi)))) +
  xlab("Sample size") +
  labs(color = "N periods", fill = "N periods") +
  theme(axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 14))+
  facet_wrap(~n_periods)

dev.off()

jpeg("plri_mse_fig.jpeg",
     width = 2000,
     height = 1500,
     res =  300)


# Line graph for mse
ggplot(plotter, aes(x = n_ind, y = sqrt(mse),
                    color = factor(n_periods), fill = factor(n_periods)
)) +
  # geom_hline(yintercept = 0, lty=2) +
  # scale_y_continuous(limits = c(0, .1)) +
  geom_line() +
  # geom_ribbon(
  #   aes(xmax = n_ind, ymin = precision_lwr, ymax = precision_upr,
  #       color = NULL),
  #   alpha = 0.15) +
  ylab("MSE") +
  xlab("Sample size")+
  labs(color = "N periods", fill = "N periods")

dev.off()


