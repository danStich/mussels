# Library load ----
library(R2jags)
library(tidyverse)

# File read ----
# List the files
files <- c("phi_dot_p_dot_50k.rda", 
           "phi_dot_p_j_50k.rda", 
           "phi_dot_p_s_50k.rda", 
           "phi_dot_p_t_50k.rda", 
           "phi_dot_p_tj_50k.rda", 
           "phi_dot_p_ts_50k.rda", 
           "phi_dot_p_tsj_50k.rda", 
           "phi_j_p_s_50k.rda", 
           "phi_js_p_s_50k.rda", 
           "phi_s_p_s_50k.rda", 
           "phi_t_p_s_50k.rda", 
           "phi_t_p_s_50k.rda", 
           "phi_tj_p_s_50k.rda", 
           "phi_tjs_p_s_50k.rda", 
           "phi_ts_p_s_50k.rda")

# Read all the files in one at a time
for(i in 1:length(files)){
  load(paste0("sabine/results/", files[i]))
}

# Detection model selection ----
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

# Survival model selection ----
s_mods <- data.frame(
  Model = c("Phi_t_p_s", "Phi_j_p_s", "Phi_s_p_s", "Phi_tj_p_s", 
            "Phi_ts_p_s", "Phi_js_p_s", "Phi_tsj_p_s"),
  DIC = c(phi_t_p_s$BUGSoutput$DIC,
          phi_j_p_s$BUGSoutput$DIC,
          phi_s_p_s$BUGSoutput$DIC,
          phi_tj_p_s$BUGSoutput$DIC,
          phi_ts_p_s$BUGSoutput$DIC,
          phi_js_p_s$BUGSoutput$DIC,
          phi_tjs_p_s$BUGSoutput$DIC))

s_mods[order(s_mods$DIC), ]


# Model selection table ----
dic_table <- rbind(det_mods[order(det_mods$DIC), ], 
                   s_mods[order(s_mods$DIC), ])

write.table(dic_table, "sabine/results/sabine-model-selection-table.csv",
            row.names = FALSE, quote = FALSE, sep = ",")

