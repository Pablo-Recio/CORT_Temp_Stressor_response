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
  mutate(trt = factor(trt,
                        levels = c("A_23", "B_23", "A_28", "B_28"),
                        labels = c("A_23" = "Control-Cold",
                                  "B_23" = "CORT-Cold",
                                  "A_28" = "Control-Hot",
                                  "B_28" = "CORT-Hot"))) %>%
  dplyr::select(-notes) %>%
data.frame()

####################################
# B) Merging mass.csv with clean individual data (clean_ind)
mass <- read.csv(here("data", "mass.csv")) %>%
  mutate(lizard_id = as.character(lizard_id))%>%
data.frame() #Read in mass data and transform lizard_id to character
mass_clean <- merge(clean_ind, mass, by = "lizard_id", all = TRUE) %>%
  mutate(delta_mass_rescaled = (delta_mass - min(delta_mass, na.rm = TRUE) + 1)) %>%
  mutate(food_ingested = food_ingested - mean(food_ingested, na.rm = TRUE))
data.frame() #Merge mass data with clean_ind data and rescale delta_mass

write.csv(mass_clean, here("./output/data_clean/mass_clean.csv")) 

####################################
# C) Merging behaviour.csv with clean individual data (clean_ind)
behav <- read.csv(here("data", "behaviour.csv"))%>%
  mutate(lizard_id = as.character(lizard_id))%>%
  mutate(day = day - 1) %>% # Adjust day to start from 0
data.frame() #Read in behaviour data and transform lizard_id to character
behav_clean <- merge(clean_ind, behav, by = "lizard_id", all = TRUE) 

write.csv(behav_clean, here("./output/data_clean/behav_clean.csv"))
