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
#
#
#| label: dataload
#| echo: false
#| warning: false
 
# Load data
	data  <-  read.csv("./data/Stressor.csv")
#
#
#
#| label: dataexplore
#| echo: false
#| warning: false

# Divide Treatment in CORT and Temperature Transform mass into mg
data <- data %>%
      mutate(temp = gsub("[AB]_", "", trt),
          cort = gsub("_[2][38]", "", trt))  %>% 
      mutate(mass=mass*1000) %>%
data.frame()
data$cort <- factor(data$cort, 
        levels = c("A" ,"B"),
        labels = c("Control", "CORT"))
data$temp <- factor(data$temp,
        levels = c("23", "28"),
        labels = c("Cold", "Hot"))
data$Day <- factor(data$Day,
        levels = c("Before", "After"))  
     
# MASS
# Create new database to analyse increases in mass (i.e. having mass after and before in different columns).
data_2 <- pivot_wider(data, id_cols=c(lizard,temp,cort,species,food_ingested),
 names_from = Day, values_from = mass)

# Estimate mass increase and rescale it to avoid negative values
data_2 <- data_2 %>% 
 mutate(delta_mass=After-Before) %>%
 mutate(rescaled_delta_mass = ifelse(is.na(delta_mass),
  NA, 
  delta_mass + abs(min(delta_mass, na.rm = T)) +1 ))
write.csv(data_2, file = "./data/Stressor_2.csv")

# Split new database by species
deli_stressor_2 <- data_2 %>% filter (species == "delicata")
guich_stressor_2 <- data_2 %>% filter (species == "guichenoti")

# Explore distribution of mass increase
  # A) L.deli
plot_mass <- ggplot(deli_stressor_2, aes(x = rescaled_delta_mass, fill = interaction(as.factor(temp),as.factor(cort), sep="_"))) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30, color = "black") +
  scale_fill_manual(values = c("28_A" = "red", "23_A" = "blue", "28_B" = "salmon", "23_B" = "lightblue")) +
  labs(x = "Mass", y = "Density", fill = "temp ~ cort") +
  theme_minimal() +
  facet_grid(as.factor(deli_stressor_2$cort)~as.factor(deli_stressor_2$temp), scales = "free_y")
print(plot_mass)

mass_normal_test<-by(deli_stressor_2$rescaled_delta_mass, interaction (deli_stressor_2$cort, deli_stressor_2$temp), function(x) shapiro.test((x)))
print(mass_normal_test)

  #B) L.guich
plot_mass <- ggplot(guich_stressor_2, aes(x = rescaled_delta_mass, fill = interaction(as.factor(guich_stressor_2$temp),as.factor(guich_stressor_2$cort), sep="_"))) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30, color = "black") +
  scale_fill_manual(values = c("28_A" = "red", "23_A" = "blue", "28_B" = "salmon", "23_B" = "lightblue")) +
  labs(x = "Mass", y = "Density", fill = "Temperature ~ CORT") +
  theme_minimal() +
  facet_grid(guich_stressor_2$cort~guich_stressor_2$temp, scales = "free_y")
print(plot_mass)

mass_normal_test<-by(guich_stressor_2$rescaled_delta_mass, interaction (guich_stressor_2$cort, guich_stressor_2$temp), function(x) shapiro.test((x)))
print(mass_normal_test)



# BEHAVIOUR
# Split data by species in original database
deli_stressor <- data %>% filter (species == "delicata")
guich_stressor <- data %>% filter (species == "guichenoti")
# Explore distribution of latency to move/time frozen
  # A) L. deli
plot_freeze <- ggplot(deli_stressor, aes(x = log(lat_to_move+1), fill = interaction(as.factor(temp),as.factor(cort), sep="_"))) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30, color = "black") +
  scale_fill_manual(values = c("28_A" = "red", "23_A" = "blue", "28_B" = "salmon", "23_B" = "lightblue")) +
  labs(x = "Mass", y = "Density", fill = "Temperature ~ CORT") +
  theme_minimal() +
  facet_grid(cort*temp~Day, scales = "free_y")
print(plot_freeze)

freeze_normal_test<-by(deli_stressor$lat_to_move, interaction (deli_stressor$cort, deli_stressor$temp, deli_stressor$Day), function(x) shapiro.test(log(x+1)))
print(freeze_normal_test)

 #B) L.guich
plot_freeze <- ggplot(guich_stressor, aes(x = log(lat_to_move+1), fill = interaction(as.factor(temp),as.factor(cort), sep="_"))) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30, color = "black") +
  scale_fill_manual(values = c("28_A" = "red", "23_A" = "blue", "28_B" = "salmon", "23_B" = "lightblue")) +
  labs(x = "Mass", y = "Density", fill = "Temperature ~ CORT") +
  theme_minimal() +
  facet_grid(cort*temp~Day, scales = "free_y")
print(plot_freeze)

freeze_normal_test<-by(guich_stressor$lat_to_move, interaction (guich_stressor$cort, guich_stressor$temp, guich_stressor$Day), function(x) shapiro.test(log(x+1)))
print(freeze_normal_test)

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
