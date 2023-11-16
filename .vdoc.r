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

pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, MASS)

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

# Divide Treatment in CORT and Temperature
data <- data %>%
      mutate(temp = gsub("[AB]_", "", trt),
          cort = gsub("_[2][38]", "", trt))  %>% 
  data.frame()

# MASS
# Create new database to analyse increases in mass (i.e. having mass after and before in different columns)
data_2 <- pivot_wider(data, id_cols=lizard_id, names_from = Day, values_from = mass)
write.csv(data_2, "./data/Stressor_2.csv")
# Transform mass into mg in new database
data_2$mass <- data$mass*1000

# Estimate mass increase

# Split data by species in new database
deli_stressor_2 <- data_2 %>% filter (species == "delicata")
guich_stressor_2 <- data_2 %>% filter (species == "guichenoti")

# Explore distribution of mass increase
  # A) L.deli
plot_mass <- ggplot(deli_stressor_2, aes(x = mass, fill = interaction(as.factor(temp),as.factor(cort), sep="_"))) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30, color = "black") +
  scale_fill_manual(values = c("28_A" = "red", "23_A" = "blue", "28_B" = "salmon", "23_B" = "lightblue")) +
  labs(x = "Mass", y = "Density", fill = "Temperature ~ CORT") +
  theme_minimal() +
  facet_grid(cort*temp~Day, scales = "free_y")
print(plot_mass)

mass_normal_test<-by(deli_stressor_2$mass, interaction (deli_stressor$cort, deli_stressor$temp, deli_stressor$Day), function(x) shapiro.test((x)))
print(mass_normal_test)

  #B) L.guich
plot_mass <- ggplot(guich_stressor_2, aes(x = mass, fill = interaction(as.factor(temp),as.factor(cort), sep="_"))) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30, color = "black") +
  scale_fill_manual(values = c("28_A" = "red", "23_A" = "blue", "28_B" = "salmon", "23_B" = "lightblue")) +
  labs(x = "Mass", y = "Density", fill = "Temperature ~ CORT") +
  theme_minimal() +
  facet_grid(cort*temp~Day, scales = "free_y")
print(plot_mass)

mass_normal_test<-by(guich_stressor$mass, interaction (guich_stressor$cort, guich_stressor$temp, guich_stressor$Day), function(x) shapiro.test((x)))
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
#| label: models1
#| echo: false
#| warning: false

# Generalized Mixed Models delicata
GMM_deli <- glmer(FC_associative ~ Associative_Trial*temp*cort + group*Associative_Trial + (1 + Associative_Trial|lizard_id), data = deli_associative, family = binomial(link = "logit"))
summary(GMM_deli)

emmeans(GMM_deli, pairwise ~ temp*Associative_Trial)
emmeans(GMM_deli, pairwise ~ cort*Associative_Trial)

# Generalized Mixed Models guichenoti
GMM_guich <- glmer(FC_associative ~ Associative_Trial*temp*cort + group*Associative_Trial + (1 + Associative_Trial|lizard_id), data = guich_associative, family = binomial(link = "logit"))
summary(GMM_guich)

emmeans(GMM_guich, pairwise ~ temp*Associative_Trial)
emmeans(GMM_guich, pairwise ~ cort*Associative_Trial)
emmeans(GMM_guich, pairwise ~ temp*cort*Associative_Trial)

#
#
#
#| label: models2
#| echo: false
#| warning: false

# Bayesian model delicata
  deli_mod <- brm(FC_associative ~ cort*temp*Associative_Trial + group*Associative_Trial + (1 + Associative_Trial|lizard_id), data = deli_associative, family = bernoulli(link = "logit"), chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
  summary(deli_mod)
  plot(deli_mod)

# Bayesian model guichenoti
  guich_mod <- brm(FC_associative ~ cort*temp*Associative_Trial + group*Associative_Trial + (1 + Associative_Trial|lizard_id), data = guich_associative, family = bernoulli(link = "logit"), chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
  summary(deli_mod)
# R2
bayes_R2(mod)

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
#| label: plots
#| echo: false
#| warning: false
 
# Plot delicata
plot_deli <- ggplot(deli_associative, aes(x = Associative_Trial, y = FC_associative, color=interaction(cort,temp))) +  
  geom_smooth(method = "glm", formula = y ~ x, method.args = list(family = "binomial"),na.action = na.exclude) +
  facet_grid(as.factor(deli_associative$group)) + 
  labs(y = "Colour Association Task", x = "Trial")  
print(plot_deli)

# Plot guichenoti
plot_guich <- ggplot(guich_associative, aes(x = Associative_Trial, y = FC_associative, color=interaction(cort,temp))) +  
  geom_smooth(method = "glm", formula = y ~ x, method.args = list(family = "binomial"),na.action = na.exclude) +
  facet_grid(as.factor(guich_associative$group)) + 
  labs(y = "Colour Association Task", x = "Trial") 
print(plot_guich)
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
