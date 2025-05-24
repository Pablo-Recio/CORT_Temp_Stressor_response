####################################
# Exploratory analyses
####################################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes,
ggh4x, cowplot, fitdistrplus, MASS, goftest)
#
source(here("R", "func.R"))

####################################
# A) Exploring distribution of the main variables for mass_clean 
data_expl <- mass_clean %>%
       mutate(treatment = paste(cort, temp, sep = "-"))
# Histogram
bin_width <- (max(data_expl$delta_mass_rescaled, na.rm = TRUE) - min(data_expl$delta_mass_rescaled, na.rm = TRUE)) / sqrt(length(data_expl$delta_mass_rescaled))
hist_plot_mass <- ggplot(data_expl, aes(x = delta_mass_rescaled)) +
  geom_histogram(binwidth = bin_width, fill = "#062d00", color = "black", alpha = 0.5) +
  theme_classic() +
  labs(x = "delta_mass", y = "Counts") +
  theme(axis.title = element_text(size = 12, family = "Times"),
        axis.text = element_text(size = 10, family = "Times"))
ggsave(here("./output/figures/exploratory/hist_mass.png"), hist_plot_mass, width = 20, height = 20, dpi = 300)

# Distribution
qq_plot_mass <- qq_plots_single(data_expl, "delta_mass_rescaled", "mass")
ggsave(here("./output/figures/exploratory/qq_plot_mass.png"), qq_plot_mass, width = 20, height = 20, dpi = 300)
#### Delta mass considered as log-normal distribution

####################################
# B) Exploring distribution of the main variables for behav_clean
