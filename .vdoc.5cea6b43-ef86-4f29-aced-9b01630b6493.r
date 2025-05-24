#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: setup

pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, MASS, dplyr, broom)

source(here("R", "func.R"))
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: models
#| echo: false
#| warning: false

# Delta mass L. delicata
  deli_mod <- brm(rescaled_delta_mass ~ food_ingested+cort*temp, data = deli_stressor_2, family = gaussian(link = "identity"), chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
  plot(deli_mod)
  summary(deli_mod)
# Delta mass L. guichenoti
  guich_mod <- brm(rescaled_delta_mass ~ food_ingested+cort*temp, data = guich_stressor_2, family = gaussian(link = "identity"), chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
  summary(guich_mod)
bayes_R2(guich_mod)

      # Seems that only food ingestion is significant in both species


# Frozen behaviour L. delicata
  deli_mod <- brm(log(lat_to_move+1) ~ Day*cort*temp + (1|lizard), data = deli_stressor, family = gaussian(link = "identity"), chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
  summary(deli_mod)
  plot(deli_mod)
bayes_R2(deli_mod)

# Frozen behaviour L. delicata
 guich_mod <- brm(log(lat_to_move+1) ~ Day*cort*temp + (1|lizard), data = guich_stressor, family = gaussian(link = "identity"), chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
  summary(guich_mod)
  plot(guich_mod)
bayes_R2(guich_mod)


# Extract posteriors
posterior <- posterior_samples(mod, pars = "^b")
plyr::ldply(lapply(posterior, mean))
plyr::ldply(lapply(posterior, function(x) quantile(x, c(0.025, 0.975))))
temp23 <- c(posterior[,"b_Intercept"], (posterior[,"b_Intercept"]+ posterior[,"b_trtB_23"]))
temp28 <- c((posterior[,"b_Intercept"] + posterior[,"b_trtB_28"]), posterior[,"b_trtB_23"])

pmcmc(temp23 - temp28)
dharma(GMM_deli2)


#
#
#
#| label: plotsMass
#| echo: false
#| warning: false
 
# Plot delta Mass L. delicata
plot_deli <- ggplot(deli_stressor_2, aes(x = food_ingested, y = rescaled_delta_mass, color=interaction(cort,temp))) +  
  stat_smooth(method = "glm", formula = y ~ x, se = TRUE, linewidth=2, alpha = 0.1) +
  scale_color_manual(values = c("CORT.Cold"="black", "Control.Cold"="grey", "CORT.Hot"="orange2", "Control.Hot"="moccasin"),
    labels=c("CORT-Cold (n=12)", "Control-Cold (n=12)", "CORT-Hot (n=11)", "Control-Hot (n=12)")
  ) +
  labs(y = "Increase in mass (mg)", x = "Ingestion rate (n items)", color="Treatment") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18, family = "sans"),  
    axis.text = element_text(size = 10, family = "sans"),  
    legend.title = element_text(size = 14, family = "sans"), 
    legend.text = element_text(size = 12, family = "sans")  
  )
print(plot_deli)
ggsave("./output/figures/fig_deli_Mass.png", plot=plot_deli, width = 18, height = 9, units = "cm", dpi = 3000)
# Plot delta Mass L. guichenoti
plot_guich <- ggplot(guich_stressor_2, aes(x = food_ingested, y = rescaled_delta_mass, color=interaction(cort,temp))) +  
  stat_smooth(method = "glm", formula = y ~ x, se = TRUE, linewidth=2, alpha = 0.1) +
  scale_color_manual(values = c("CORT.Cold"="black", "Control.Cold"="grey", "CORT.Hot"="orange2", "Control.Hot"="moccasin"),
    labels=c("CORT-Cold (n=10)", "Control-Cold (n=7)", "CORT-Hot (n=11)", "Control-Hot (n=10)")
  ) +
  labs(y = "Increase in mass (mg)", x = "Ingestion rate (n items)", color="Treatment") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18, family = "sans"),  
    axis.text = element_text(size = 10, family = "sans"),  
    legend.title = element_text(size = 14, family = "sans"), 
    legend.text = element_text(size = 12, family = "sans")  
  )
print(plot_guich)
ggsave("./output/figures/fig_guich_Mass.png", plot=plot_guich, width = 18, height = 9, units = "cm", dpi = 3000)
#
#
#
#| label: plotsFrozenbehaviour
#| echo: false
#| warning: false

# Plot time frozen L. delicata
plot_deli_BEH <- ggplot(deli_stressor, aes(x = Day, y = log(lat_to_move+1), fill=interaction(cort,temp))) +  
  geom_boxplot(outlier.shape = NA) +
  ylim(0, 10) +
  scale_fill_manual(values = c("CORT.Cold"="black", "Control.Cold"="grey", "CORT.Hot"="orange2", "Control.Hot"="moccasin"),
    labels=c("CORT-Cold (n=12)", "Control-Cold (n=12)", "CORT-Hot (n=11)", "Control-Hot (n=12)")
  ) +
  labs(y = "Time frozen (log-transformed)", x = "Day", fill="Treatment") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18, family = "sans"),  
    axis.text = element_text(size = 12, family = "sans"),  
    legend.title = element_text(size = 14, family = "sans"), 
    legend.text = element_text(size = 12, family = "sans")  
  )
print(plot_deli_BEH)
ggsave("./output/figures/fig_deli_BEH.png", plot=plot_deli_BEH, width = 15, height = 12, units = "cm", dpi = 3000)

# Plot time frozen L. guichenoti
plot_guich_BEH <- ggplot(guich_stressor, aes(x = Day, y = log(lat_to_move+1), fill=interaction(cort,temp))) +  
  geom_boxplot(outlier.shape = NA) +
  ylim(0, 8) +
  scale_fill_manual(values = c("CORT.Cold"="black", "Control.Cold"="grey", "CORT.Hot"="orange2", "Control.Hot"="moccasin"),
    labels=c("CORT-Cold (n=10)", "Control-Cold (n=7)", "CORT-Hot (n=11)", "Control-Hot (n=10)")
  ) +
  labs(y = "Time frozen (s)", x = "Day", fill="Treatment") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18, family = "sans"),  
    axis.text = element_text(size = 12, family = "sans"),  
    legend.title = element_text(size = 14, family = "sans"), 
    legend.text = element_text(size = 12, family = "sans")  
  )
print(plot_guich_BEH)
ggsave("./output/figures/fig_guich_BEH.png", plot=plot_guich_BEH, width = 15, height = 12, units = "cm", dpi = 3000)

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: fig-fig1
#| fig-cap: "Activity budget of the focal group of chacma baboons (Papio ursinus) at the Cape Peninsula, South Africa, during the study period. The data are presented as the mean percentage of time spent in each activity category per hour (± SE)."

data %>% ggplot(aes(x = Lizard_id))  + geom_histogram()


#
#
#
#
#
#| label: tbl-tb1
#| tbl-cap: "Summary of the data"

tab <- data %>% 
  group_by(Lizard_id) %>% 
  summarise(
    n = n()
  ) 
flextable(tab[2:3,])
#
#
#
#
#
