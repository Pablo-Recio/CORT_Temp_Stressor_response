#################### 
####################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom,
tidybayes, ggh4x, cowplot, fitdistrplus, MASS, goftest)
#
####################
####################
# Extract sample size per group
#' @title sample
#' @description Estract sample size per group
#' @param df Database
#' @param corti To select cort treatment ("CORT"/"Control")
#' @param therm To select temp treatment ("Cold"/"Hot")
sample <- function(df, corti, therm){
  sample_size <- df %>%
    filter(cort == corti, temp == therm) %>%
    distinct(lizard_id) %>%
    summarise(total_count = n()) %>%
    pull(total_count)
  return(sample_size)
}
####################
####################
# Function to generate four Q-Q plots for one variable
#' @title qq_plots_single 
#' @description Generate four Q-Q plots for one variable
#' @param df Database
#' @param vle Variable chosen
#' @param label Label assigned to each variable
#' @return fig_qq
qq_plots_single <- function(df, vle, label) {
  # Remove NAs
  plot_df <- data.frame(value = na.omit(df[[vle]])) 

    # Normal Q-Q Plot
    p1 <- ggplot(plot_df, aes(sample = value)) +
      stat_qq() + 
      stat_qq_line(color = "red") +
      ggtitle(paste(label, "- Normal Q-Q Plot")) +
      theme_classic()

    # Log-Normal Q-Q Plot
    p2 <- ggplot(plot_df, aes(sample = log(value +1))) +
      stat_qq() + 
      stat_qq_line(color = "red") +
      ggtitle(paste(label, "- Log-Normal Q-Q Plot")) +
      theme_classic()

    # Exponential Q-Q Plot
    p3 <- ggplot(plot_df, aes(sample = value)) +
      stat_qq(distribution = qexp, dparams = list(rate = 1/mean(plot_df$value))) +
      stat_qq_line(distribution = qexp, dparams = list(rate = 1/mean(plot_df$value)), color = "red") +
      ggtitle(paste(label, "- Exponential Q-Q Plot")) +
      theme_classic()

    # Arrange plots in a grid
    fig_qq <- plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(0.33, 0.34, 0.33))
    
  return(fig_qq)
}
####################
####################
# Fit models and extract posteriors both per each treatment:
#' @title fit_m Function
#' @description Fit brm models for different the associative task
#' @param df To select the df
#' @param cat To label whether the models are the preliminary or not
#' @param var To select the variable
#' @param formula To select the most appropriate formula
#' @param fam To select the most appropriate distribution depending on the main variable
#' @param label To select the most appropriate label (OB or OT)
#' @param refit To choose whether to refit the models (TRUE, default) or use the ones already made (FALSE)
#' @return Raw posteriors of fitted brm model for each treatment, species, and group (df)
fit_m <- function(df, cat, var, formula, fam, label, refit = TRUE) {
  #Fit the model only if it has not been fit yet (if refit=TRUE)
  if(refit){
    # Fit the model
    model <- brm(formula,
                data = df,
                family = fam,
                chains = 4, cores = 4, iter = 8000, warmup = 2000,
                control = list(adapt_delta = 0.99, max_treedepth = 11))
    # Write the model to a file
    saveRDS(model, file = paste0(here("output/models/"), var, "_", cat, "_", label, ".rds"))
  } else {
      # Read the model from a file
      model <- readRDS(file = paste0(here("output/models/"), var, "_", cat, "_", label, ".rds"))
  } 
  # Extract posteriors
  posteriors <- as_draws_df(model)
  return(posteriors)
}
####################
####################
# Estimate pmcm
#' @title pMCMC Function
#' @param x The vector for the posterior distribution. Note that this will test the null hypothesis that the parameter of interest is significantly different from 0. 
#' @param null A numeric value decsribing what the null hypothesis should be
#' @param twotail Whether to conduct a one-tailed hypothesis or a two-tailed hypotheses. Default = true indicating a two-tailed test will be done.
#' @param dir The direction of the one-tail test (<) or (>)
pmcmc <- function(x, null = 0, twotail = TRUE, dir){
  if(twotail){
    1.645*(1 - max(table(x<=null) / length(x)))
  } else{
    if(dir == "<"){
      (max(sum(x>=null) / length(x)))
    } else{
      if(dir == ">"){
        (max(sum(x<=null) / length(x)))
      } else{
        stop("dir not valid")
      }
    }
  }
}
####################
####################
# Function to format numbers with n decimal places
#' @title format_dec
#' @param x The object
#' @param n The number of decimals
format_dec <- function(x, n) {
  z <- sprintf(paste0("%.",n,"f"), x)
  return(z)
}
####################
####################
# Function to format p_values with n decimal places
#' @title format_p
#' @param x The object
#' @param n The number of decimals
#' @param equal To write whether we want the equal or not
format_p <- function(x, n, equal) {
  if (equal == TRUE){
    label <- "="
  } else {
    label <- ""
  }
  z <- sprintf(paste0("%.",n,"f"), x)
  tmp <- ifelse(as.numeric(z) <= 0.001, "< 0.001",
         ifelse(as.numeric(z) <= 0.05 & as.numeric(z) > 0.001, "< 0.05",
                paste0(label, as.character(z))))
  return(tmp)
}
###################
###################
# Function to make the first cleaning for the posteriors and create a big df that includes
# all variables and their predictors
#' @title tidy_post
#' @param df to select the posteriors to use. Options: choice or errors
tidy_post <- function(df) {
  df <- as.data.frame(df)
  summary_df <- df %>%
    dplyr::select(starts_with("b_")) %>%  # Select only parameter estimates
    pivot_longer(everything(), names_to = "Predictor", values_to = "Estimate") %>%
    group_by(Predictor) %>%
    summarise(
      Mean = format_dec(mean(Estimate), 3),
      CI_lower = format_dec(quantile(Estimate, 0.05, na.rm = TRUE), 3),
      CI_upper = format_dec(quantile(Estimate, 0.95, na.rm = TRUE), 3),
      PMCMC = format_dec(pmcmc(Estimate, null = 0, twotail = TRUE), 3),
      .groups = "drop"
    )
  
  return(summary_df)
}
####################
####################
# Function to rename and refine the posterior databases for preparing the tables.
#' @title refine_post
#' @param df to select the databse
refine_post <- function(df) {
  variable_order <- c("Mit density",
                     "Mit potential",
                     "ROS",
                     "DNA damage",
                     "Peroxidation",
                     "Detection lat")
  predictor_order <- c("b_Intercept",
                       "b_cortCORT",
                       "b_tempHot",
                       "b_cortCORT:tempHot",
                       "b_age_euthanasia",
                       "b_sexFemale",
                       "b_motivation",
                       "b_motivation:cortCORT")
  df <- df %>%
    mutate(
      Model = case_when(
        Model == "m_def_mean_mitodensity_OB" ~ "Mit density",
        Model == "m_def_mean_potential_OB" ~ "Mit potential",
        Model == "m_def_mean_ros_OB" ~ "ROS",
        Model == "m_def_mean_dnadamage_OB" ~ "DNA damage",
        Model == "m_def_mean_peroxidation_OB" ~ "Peroxidation",
        Model == "m_def_t_D_Chemical" ~ "Detection lat",
        Model == "m_def_mean_mitodensity_OT" ~ "Mit density",
        Model == "m_def_mean_potential_OT" ~ "Mit potential",
        Model == "m_def_mean_ros_OT" ~ "ROS",
        Model == "m_def_mean_dnadamage_OT" ~ "DNA damage",
        Model == "m_def_mean_peroxidation_OT" ~ "Peroxidation",
        Model == "m_def_t_D_Visual" ~ "Detection lat",
        TRUE ~ Model
      )
    ) %>%
    mutate(CI_95 = paste0("[", CI_lower, " , ", CI_upper, "]")) %>%
    rename(
      Variable = Model,
      Predictors = Predictor,
      'Estimate Mean' = Mean, 
      '95% CI' = CI_95, 
      PMCMC = PMCMC
    ) %>%
    dplyr::select(
      Variable,
      Predictors,
      'Estimate Mean',
      '95% CI',
      PMCMC) %>%
    mutate(Variable = factor(Variable, levels = variable_order)) %>%
    mutate(Predictors = factor(Predictors, levels = predictor_order))
  
  return(df)
}
####################
####################
# Function to get estimates to make the contrasts between treatments
#' @title post_values
#' @param df To select the df
#' @param fac To indicate if there are other parameters to take into account
post_values <- function(df, fac){
  df <- as.data.frame(df)
  #Getting the estimated values per each treatment
  if(fac == "none"){
    Control_Cold <- df$b_Intercept
    Control_Hot <- df$b_Intercept + df$b_tempHot
    CORT_Cold <- df$b_Intercept + df$b_cortCORT
    CORT_Hot <- df$b_Intercept + df$b_cortCORT + df$b_tempHot + df$`b_cortCORT:tempHot`
  } else if (fac == "sex"){
    # Males
    Control_Cold_males <- base::sample(df$b_Intercept, size = 12000, replace = FALSE)
    Control_Hot_males <- base::sample(df$b_Intercept + df$b_tempHot, size = 12000, replace = FALSE)
    CORT_Cold_males <- base::sample(df$b_Intercept + df$b_cortCORT, size = 12000, replace = FALSE)
    CORT_Hot_males <- base::sample(df$b_Intercept + df$b_cortCORT + df$b_tempHot +
                                      df$`b_cortCORT:tempHot`, size = 12000, replace = FALSE)
    # Females
    Control_Cold_fem <- base::sample(df$b_Intercept + df$b_sexFemale, size = 12000, replace = FALSE)
    Control_Hot_fem <- base::sample(df$b_Intercept + df$b_tempHot + df$b_sexFemale, size = 12000, replace = FALSE)
    CORT_Cold_fem <- base::sample(df$b_Intercept + df$b_cortCORT + df$b_sexFemale, size = 12000,
                                    replace = FALSE)
    CORT_Hot_fem <- base::sample(df$b_Intercept + df$b_cortCORT + df$b_tempHot +
                                    df$`b_cortCORT:tempHot` + df$b_sexFemale, size = 12000, replace = FALSE)
    # Pool
    CORT_Cold <- c(CORT_Cold_males, CORT_Cold_fem)
    CORT_Hot <- c(CORT_Hot_males, CORT_Hot_fem)
    Control_Cold <- c(Control_Cold_males, Control_Cold_fem)
    Control_Hot <- c(Control_Hot_males, Control_Hot_fem)
  }
  data_values <- data.frame(
    CORT_Cold = CORT_Cold,
    CORT_Hot = CORT_Hot,
    Control_Cold = Control_Cold,
    Control_Hot = Control_Hot
  )
  return(data_values)
}
####################
####################
# Function to create the plots
#' @title plotting
#' @param df to select the df (only posteriors)
#' @param lab to select the appropriate label for axix Y
plotting <- function(df, lab){
  # Build db for violin plots
  db_violin <- data.frame(
    Treatment = rep(c("CORT-Cold (n = 20)",
                      "CORT-Hot (n = 20)",
                      "Control-Cold (n = 20)",
                      "Control-Hot (n = 20)"), each = nrow(df)),
    Estimate = c(df$CORT_Cold, df$CORT_Hot, df$Control_Cold, df$Control_Hot)
  ) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Control-Hot (n = 20)",
                                       "CORT-Hot (n = 20)",
                                       "Control-Cold (n = 20)",
                                       "CORT-Cold (n = 20)")))

  # Build db for geom_bars
  db_bars <- db_violin %>%
    group_by(Treatment) %>%
    summarize(
      Mean = mean(Estimate),
      SD = sd(Estimate),
      SE = sd(Estimate)/sqrt(length(Estimate))
    )
  # Making the plot
  plot <- ggplot(db_violin, aes(x = Treatment, y = Estimate, fill = Treatment)) +
  geom_flat_violin(alpha = 0.5) +
  scale_fill_manual(values = c("Control-Hot (n = 20)" = "#fa927d",
                              "CORT-Hot (n = 20)" = "#b50101",
                              "Control-Cold (n = 20)" = "#68bde1",
                              "CORT-Cold (n = 20)" = "#00008B")) +
  geom_point(data = db_bars, aes(y = Mean, x = Treatment), position = position_dodge(width = 0.75), color = "black", fill = "black", size = 3) +
  geom_segment(data = db_bars, aes(y = Mean - SD, yend = Mean + SD, x = Treatment, xend = Treatment), size = 1.5, color = "black") +
  ylim(min(db_violin$Estimate), max(db_violin$Estimate)) +
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic() +
  labs(y = lab, x = "Treatments") +
  theme(plot.margin = margin(3, 3, 3, 3, "mm"),
    axis.text.y = element_blank(),   # Remove y-axis labels
    axis.title = element_text(size = 12, family = "Times"),
    axis.text = element_text(size = 10, family = "Times"),
    legend.position = "none",
    legend.title = element_text(size = 12, family = "Times"),
    legend.text = element_text(size = 11, family = "Times")
    )
  return(plot)
}
