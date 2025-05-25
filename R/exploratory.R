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
data_expl <- behav_clean %>%
  mutate(treatment = paste(cort, temp, sep = "-")) %>%
  mutate(
    lat_move = as.numeric(lat_move),
    lat_shelter = as.numeric(lat_shelter)
  )

# Variables to explore
var <- c("lat_move", "lat_shelter")

# Histograms
hist_list <- list()
for (i in var) {
  bin_width <- (max(data_expl[[i]], na.rm = TRUE) - min(data_expl[[i]], na.rm = TRUE)) / sqrt(length(data_expl[[i]]))
  # Labels
  if (i == "lat_move"){
     lab <- "Latency to move (s)"
  } else if (i == "lat_shelter") {
     lab <- "Latency to shelter (s)"
  }
  # Create histogram
  hist_plot <- ggplot(data_expl, aes(x = .data[[i]])) +
    geom_histogram(binwidth = bin_width, fill = "#062d00", color = "black", alpha = 0.5) +
    theme_classic() +
    labs(x = lab, y = "Counts") +
    theme(axis.title = element_text(size = 12, family = "Times"),
          axis.text = element_text(size = 10, family = "Times"))
  hist_list[[i]] <- hist_plot
}
hist_behav <- plot_grid(plotlist = hist_list, ncol = 3)
ggsave(here("./output/figures/exploratory/hist_behav.png"), hist_behav, width = 20, height = 20, dpi = 300)

# Distribution
qq_plots_list <- lapply(var, function(v) {
  if (i == "lat_move"){
     label <- "Latency to move (s)"
  } else if (i == "lat_shelter") {
     label <- "Latency to shelter (s)"
  } else {
    v  # Fallback to variable name if not found
  }
  qq_plots_single(data_expl, v, label)
})
#
behav_qq_plot <- plot_grid(plotlist = qq_plots_list, ncol = 1)
ggsave(here("./output/figures/exploratory/behav_qq_plot.png"), behav_qq_plot, width = 20, height = 20, dpi = 300)
#
#