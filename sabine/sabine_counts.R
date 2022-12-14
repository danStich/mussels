# Libraries ----
library(tidyverse)
library(lubridate)

# Data read ----
sabine <- read.csv("sabine.csv")

sabine$date <- as.Date(as.character(sabine$date), format = "%m/%d/%Y")
sabine$month <- month(sabine$date)
sabine$year <- year(sabine$date)


# Drop dupes
# inds <- sabine %>%
#   filter(!duplicated(cbind(hallprint_id)))

# Summaries ----
inds %>% 
  group_by(species, site, date) %>% 
  summarize(count = n())

inds %>% 
  group_by(species, site, month, year) %>% 
  summarize(count = n())

inds %>% 
  group_by(species, site, year) %>% 
  summarize(count = n())

# Make a dummy capture history ----
ch <- data.frame(reshape::cast(data = inds, formula = hallprint_id ~ date))
head(ch)

ch[is.na(ch)] <- 0
for(i in 2:ncol(ch)){
  ch[,i][ch[,i]>0] = 1
}

hist(rowSums(ch[,2:4]))




