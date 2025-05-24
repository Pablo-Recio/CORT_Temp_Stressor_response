####################################
# Data_process
####################################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2,
               lme4, zoo, lmerTest, broom, tidybayes, forcats)
#
source(here("R", "func.R"))
# A) Cleaning the general_data.csv for all individuals
clean_ind <- read.csv(here("data", "general_data.csv")) %>%
  mutate(lizard_id = as.character(lizard_id)) %>%
  mutate(temp = gsub("[AB]_", "", trt),
         cort = gsub("_[2][38]", "", trt))  %>%
  mutate(temp = factor(temp,
                       levels = c("23", "28"),
                       labels = c("23" = "Cold", "28" = "Hot"))) %>%
  mutate(cort = factor(cort,
                       levels = c("A", "B"),
                       labels = c("A" = "Control", "B" = "CORT"))) %>%
  dplyr::select(-notes, -trt) %>%
data.frame()

####################################
# B) Merging mass.csv with clean individual data (clean_ind)
mass <- read.csv(here("data", "mass.csv")) %>%
  mutate(lizard_id = as.character(lizard_id))%>%
data.frame() #Read in mass data and transform lizard_id to character
mass_clean <- merge(clean_ind, mass, by = "lizard_id", all = TRUE) %>%
  mutate(delta_mass_rescaled = (delta_mass - min(delta_mass, na.rm = TRUE))) %>%
data.frame() #Merge mass data with clean_ind data and rescale delta_mass

write.csv(mass_clean, here("./output/data_clean/mass_clean.csv")) 

####################################
# C) Merging behaviour.csv with clean individual data (clean_ind)
behav <- read.csv(here("data", "behaviour.csv"))%>%
  mutate(lizard_id = as.character(lizard_id))%>%
data.frame() #Read in behaviour data and transform lizard_id to character
behav_clean <- merge(clean_ind, behav, by = "lizard_id", all = TRUE) 

write.csv(behav_clean, here("./output/data_clean/behav_clean.csv"))
