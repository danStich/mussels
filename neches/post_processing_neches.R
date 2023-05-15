library(tidyverse)
library(reshape)
library(lubridate)

#setwd("~/Grad work/mussels/neches")
load("results/phi_ts_p_s_50000_3.rda")

# Data read ----
neches <- read.csv("data/neches.csv", stringsAsFactors = FALSE)

# Have a look
glimpse(neches)

# Remove individuals without a Hallprint tag
mussels <- neches %>% 
  filter(hallprint1 != "") %>% 
  mutate(individual = paste0(hallprint1, hallprint2))

# Get species for each individual
ind_spp <- data.frame(
  unique(cbind(mussels$individual[mussels$species != ""], 
               mussels$species[mussels$species != ""])))
names(ind_spp) <- c("individual", "species_2")

# Merge back with raw data to get species in "mussels"
mussels <- merge(mussels, ind_spp, by = "individual")



survivals <- melt(phi_ts_p_s$BUGSoutput$sims.list$lphi)

View(survivals)

names(survivals) <- c("iteration","time", "species", "lphi")

survivals <- survivals %>% filter(time == 1)

view(survivals)

survivals$phi <- boot::inv.logit(survivals$lphi)

unique(mussels$species)

##PLRI	Pleurobema riddellii 17
##POAM	Potamilus amphichaenus 8

# getting survival for plri & poam ----
#plir
plri <- survivals %>% 
  filter(species == 17)

#mean surivial rate for louisianna pigtoe simulation
mean(plri$phi)

# poam
poam <- survivals %>% 
  filter(species ==8)

# mean survival rate for texas heelsplitter simulation
mean(poam$phi)

# getting detection for plri & poam ----

detections <- melt(phi_ts_p_s$BUGSoutput$sims.list$lp)

view(detections)

names(detections) <- c("iteration", "species", "lphi")

head(detections)

detections$p <- boot::inv.logit(detections$lp)


plri <- detections %>% 
  filter(species == 17)

# mean detection for plri for simulations
mean(plri$p)


poam <- detections %>% 
  filter(species == 8)

# mean deteciton for poam for simulations
mean(poam$p)
