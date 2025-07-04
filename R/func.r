#################### 
####################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, ggh4x)
#
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
# Extract sample size per group
#' @title sample
#' @description Estract sample size per group
#' @param df To select the df
#' @param corti To select cort treatment ("CORT"/"Control")
#' @param therm To select temp treatment ("Cold"/"Hot")
sample <- function(df, corti, therm){
  sample_size <- df %>%
                filter(cort == corti, temp == therm) %>%
                group_by(lizard_id) %>%
                summarise(n = n()) %>%
                summarise(total_count = sum(n)) %>%
                pull(total_count)
  return(sample_size)
}
####################
####################
# Fit models and extract posteriors.
#' @title fit_m_beh Function
#' @description Fit brm models for different the associative task
#' @param df To select the df
#' @param fam To specify the family of the model
#' @param sp To specify the species
#' @param var To specify the response variable
#' @param type To specify which type of model/formula
#' @return Raw posteriors of fitted brm model for each treatment, species, and group (df)
fit_m <- function(df, fam, sp, type, var) {
  #Fit the model only if it has not been fit yet (if refit=TRUE)
  if(refit){
    if(type == "beh"){
      formula <- as.formula(paste(var, 
                            "~ day*trt + (1|lizard_id) + (1|clutch)"))
    } else if(type == "mass"){
      formula <- as.formula(paste(var, 
                            "~ trt + food_ingested + (1|clutch)"))
    }
    # Fit the model
    model <- brm(formula = formula,
                data = df,
                family = fam,
                chains = 4, cores = 4, iter = 10000, warmup = 4000,
                control = list(adapt_delta = 0.999, max_treedepth = 15))
    # Write the model to a file
    saveRDS(model, file = paste0(here("output/models/"), sp, "_", var, ".rds"))
  } else {
      # Read the model from a file
      model <- readRDS(file = paste0(here("output/models/"), sp, "_", var, ".rds"))
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
    2*(1 - max(table(x<=null) / length(x)))
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
# Function to clean the posteriors and create a big df that includes
# all variables and their predictors
#' @title tidy_post
#' @param df to select the posteriors to use. 
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
      pMCMC = format_dec(pmcmc(Estimate, null = 0, twotail = TRUE), 3),
      .groups = "drop"
    ) %>%
  data.frame()
  return(summary_df)
}
####################
####################
# Function to get estimates to make the contrasts between treatments
#' @title post_values
#' @param df To select the df
#' @param var To indicate the kind of variable to use (beh or mass)
post_values <- function(df, var){
  df <- as.data.frame(df)
  if(var == "beh"){
    int_Control_Cold <- df$b_Intercept
    int_CORT_Cold <- df$b_Intercept + df$b_trtCORTMCold
    int_Control_Hot <- df$b_Intercept + df$b_trtControlMHot
    int_CORT_Hot <- df$b_Intercept + df$b_trtCORTMHot
    slope_Control_Cold <- df$b_day
    slope_CORT_Cold <- df$b_day + df$`b_day:trtCORTMCold`
    slope_Control_Hot <- df$b_day + df$`b_day:trtControlMHot`
    slope_CORT_Hot <- df$b_day + df$`b_day:trtCORTMHot`    
    # Create a data frame with the estimates
    data_values <- data.frame(
      int_Control_Cold = int_Control_Cold,
      int_CORT_Cold = int_CORT_Cold,
      int_Control_Hot = int_Control_Hot,
      int_CORT_Hot = int_CORT_Hot,
      slope_Control_Cold = slope_Control_Cold,
      slope_CORT_Cold = slope_CORT_Cold,
      slope_Control_Hot = slope_Control_Hot,
      slope_CORT_Hot = slope_CORT_Hot)
  
  } else if (var == "mass") {
    int_Control_Cold <- df$b_Intercept
    int_CORT_Cold <- df$b_Intercept + df$b_trtCORTMCold
    int_Control_Hot <- df$b_Intercept + df$b_trtControlMHot
    int_CORT_Hot <- df$b_Intercept + df$b_trtCORTMHot
    # Create a data frame with the estimates
    data_values <- data.frame(
      int_Control_Cold = int_Control_Cold,
      int_CORT_Cold = int_CORT_Cold,
      int_Control_Hot = int_Control_Hot,
      int_CORT_Hot = int_CORT_Hot)
  }
  return(data_values)
}
####################
####################
# Function to get estimates to make the contrasts between treatments
#' @title contrasts_func
#' @param sp To select the species ("guich" or "deli")
#' @param res To select the response variable ("move", "shelter", "emergence", "mass")
contrasts_func <- function(sp, res){
  data_table <- data.frame()
  df <- get(paste0(sp, "_", res, "_posteriors"))
  if(res == "mass"){
    # Contrasts
    Temperature_int <- format_dec(mean(c(df$int_CORT_Hot, df$int_Control_Hot)) 
                                  - mean(c(df$int_CORT_Cold, df$int_Control_Cold)), 3)

    pMCMC_temp_int <- format_p(pmcmc(c(df$int_CORT_Hot, df$int_Control_Hot) 
                                    - c(df$int_CORT_Cold, df$int_Control_Cold)), 3, equal = FALSE)
    #
    CORT_int <- format_dec(mean(c(df$int_Control_Hot, df$int_Control_Cold)) 
                          - mean(c(df$int_CORT_Hot, df$int_CORT_Cold)), 3)

    pMCMC_cort_int <- format_p(pmcmc(c(df$int_Control_Hot, df$int_Control_Cold) 
                                    - c(df$int_CORT_Hot, df$int_CORT_Cold)), 3, equal = FALSE)
    #
    Interaction_int <- format_dec((mean(df$int_Control_Hot) - mean(df$int_CORT_Hot)) 
                                - (mean(df$int_Control_Cold) - mean(df$int_CORT_Cold)), 3)

    pMCMC_int_int <- format_p(pmcmc((df$int_Control_Hot - df$int_CORT_Hot) 
                              - (df$int_Control_Cold - df$int_CORT_Cold)), 3, equal = FALSE)
    data_temp <- data.frame(#Intercepts
                          mean_Temperature_int = as.numeric(Temperature_int),
                          pMCMC_Temperature_int = pMCMC_temp_int,
                          mean_Hormone_int = as.numeric(CORT_int),
                          pMCMC_Hormone_int = pMCMC_cort_int,
                          mean_Interaction_int = as.numeric(Interaction_int),
                          pMCMC_Interaction_int = pMCMC_int_int)
    data_table <- dplyr::bind_rows(data_table, data_temp)
  } else {
    # Intercept contrasts
    Temperature_int <- format_dec(mean(c(df$int_CORT_Hot, df$int_Control_Hot)) 
                                  - mean(c(df$int_CORT_Cold, df$int_Control_Cold)), 3)

    pMCMC_temp_int <- format_p(pmcmc(c(df$int_CORT_Hot, df$int_Control_Hot) 
                                    - c(df$int_CORT_Cold, df$int_Control_Cold)), 3, equal = FALSE)
    #
    CORT_int <- format_dec(mean(c(df$int_Control_Hot, df$int_Control_Cold)) 
                          - mean(c(df$int_CORT_Hot, df$int_CORT_Cold)), 3)

    pMCMC_cort_int <- format_p(pmcmc(c(df$int_Control_Hot, df$int_Control_Cold) 
                                    - c(df$int_CORT_Hot, df$int_CORT_Cold)), 3, equal = FALSE)
    #
    Interaction_int <- format_dec((mean(df$int_Control_Hot) - mean(df$int_CORT_Hot)) 
                                - (mean(df$int_Control_Cold) - mean(df$int_CORT_Cold)), 3)

    pMCMC_int_int <- format_p(pmcmc((df$int_Control_Hot - df$int_CORT_Hot) 
                              - (df$int_Control_Cold - df$int_CORT_Cold)), 3, equal = FALSE)

    # Slope contrasts
    Temperature_slope <- format_dec(mean(c(df$slope_CORT_Hot, df$slope_Control_Hot)) 
                                  - mean(c(df$slope_CORT_Cold, df$slope_Control_Cold)), 3)

    pMCMC_temp_slope <- format_p(pmcmc(c(df$slope_CORT_Hot, df$slope_Control_Hot) 
                                      - c(df$slope_CORT_Cold, df$slope_Control_Cold)), 3, equal = FALSE)
    #
    CORT_slope <- format_dec(mean(c(df$slope_Control_Hot, df$slope_Control_Cold)) 
                            - mean(c(df$slope_CORT_Hot, df$slope_CORT_Cold)), 3)

    pMCMC_cort_slope <- format_p(pmcmc(c(df$slope_Control_Hot, df$slope_Control_Cold) 
                                      - c(df$slope_CORT_Hot, df$slope_CORT_Cold)), 3, equal = FALSE)
    #
    Interaction_slope <- format_dec((mean(df$slope_Control_Hot) - mean(df$slope_CORT_Hot)) 
                                  - (mean(df$slope_Control_Cold) - mean(df$slope_CORT_Cold)), 3)

    pMCMC_int_slope <- format_p(pmcmc((df$slope_Control_Hot - df$slope_CORT_Hot) 
                                    - (df$slope_Control_Cold - df$slope_CORT_Cold)), 3, equal = FALSE)
    ##   
    data_temp <- data.frame(#Intercepts
                            mean_Temperature_int = as.numeric(Temperature_int),
                            pMCMC_Temperature_int = pMCMC_temp_int,
                            mean_Hormone_int = as.numeric(CORT_int),
                            pMCMC_Hormone_int = pMCMC_cort_int,
                            mean_Interaction_int = as.numeric(Interaction_int),
                            pMCMC_Interaction_int = pMCMC_int_int,
                            #Slopes
                            mean_Temperature_slope = as.numeric(Temperature_slope),
                            pMCMC_Temperature_slope = pMCMC_temp_slope,
                            mean_Hormone_slope = as.numeric(CORT_slope),
                            pMCMC_Hormone_slope = pMCMC_cort_slope,
                            mean_Interaction_slope = as.numeric(Interaction_slope),
                            pMCMC_Interaction_slope = pMCMC_int_slope)
    data_table <- dplyr::bind_rows(data_table, data_temp)
  }
  return(data_table)
}
####################
####################
# Function to create the plots for the behaviours
#' @title plotting_beh
#' @param sp To select the species ("guich" or "deli")
#' @param res To select the response variable ("move", "shelter", "emergence", "mass")
#' @param type To select whether the type of the plot is for the intercept ("int"),
# the slope ("slop"), or the change over trials ("trials")
plotting <- function(sp, res, type){

    # Labels for the y-axis
  if(type == "slop"){
    label <- "Estimated slopes"
  } else {   
    if(res == "move"){
      label <- "Latency to move (s)"
    } else if(res == "shelter"){
      label <- "Latency to shelter (s)"
    } else if(res == "emergence"){
      label <- "Probability of emergence"
    } else if(res == "mass"){
      label <- "Changes in mass"
    }
  }

  # Getting the df
  df <- get(paste0(sp, "_", res, "_posteriors"))
  
  # Getting sample sizes
  if(sp == "guich"){
    n_Control_Cold <- n_list_guichenoti$`Control-Cold`
    n_CORT_Cold <- n_list_guichenoti$`CORT-Cold`
    n_Control_Hot <- n_list_guichenoti$`Control-Hot`
    n_CORT_Hot <- n_list_guichenoti$`CORT-Hot`
  } else if(sp == "deli"){
    n_Control_Cold <- n_list_delicata$`Control-Cold`
    n_CORT_Cold <- n_list_delicata$`CORT-Cold`
    n_Control_Hot <- n_list_delicata$`Control-Hot`
    n_CORT_Hot <- n_list_delicata$`CORT-Hot`
  }

  # Improving the legend labels
  fill_colors <- setNames(
    c("#68bde1", "#00008B", "#b50101", "#fa927d"),
    c(paste0("Control-Cold (n=", n_Control_Cold, ")"),
      paste0("CORT-Cold (n=", n_CORT_Cold, ")"),
      paste0("CORT-Hot (n=", n_CORT_Hot, ")"),
      paste0("Control-Hot (n=", n_Control_Hot, ")")))

  # Plots for the trials
  if(type == "trials"){     
    # Create new matrix and new dfs
    treat_matrix <- matrix(NA, nrow = 24000, ncol = 8)
    fig_df <- data.frame()
    treatments <- c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot")
    for(t in treatments){
      if(t == "CORT-Cold"){
        u <- df$int_CORT_Cold
        m <- df$slope_CORT_Cold
      } else if(t == "Control-Cold"){
        u <- df$int_Control_Cold
        m <- df$slope_Control_Cold
      } else if(t == "CORT-Hot"){
        u <- df$int_CORT_Hot
        m <- df$slope_CORT_Hot
      } else if(t == "Control-Hot"){
        u <- df$int_Control_Hot
        m <- df$slope_Control_Hot
      } else {
      stop("loop wrong")
      }
      # Getting the final data for the plot based on the distribution of the response variable
      for(x in 0:7) {
        for(j in 1:24000){
          if (res == "emergence") {
            value <- plogis(u[j] + m[j] * x)
          } else {
            value <- exp(u[j] + m[j] * x)
          }
          treat_matrix[j, x + 1] <- value 
        }
      }
      treat_df <- as.data.frame(treat_matrix)
      colnames(treat_df) <- as.character(1:8)  # Adjust column names
      treat_df <- gather(treat_df, key = "Trial", value = "Value")  # Reshape data frame
      treat_df$Treatment <- t
      fig_df <- rbind(fig_df, treat_df)
    }

    plot_df <- fig_df %>%
      mutate(Treatment = factor(Treatment,
            levels = c("CORT-Hot", "Control-Hot", "CORT-Cold", "Control-Cold"),
            labels = c("CORT-Hot" = paste0("CORT-Hot (n=", n_CORT_Hot, ")"),
                      "Control-Hot" = paste0("Control-Hot (n=", n_Control_Hot, ")"),
                      "CORT-Cold" = paste0("CORT-Cold (n=", n_CORT_Cold, ")"),
                      "Control-Cold" = paste0("Control-Cold (n=", n_Control_Cold, ")")))) %>%
      mutate(Trial = as.numeric(Trial)) %>%
      group_by(Trial, Treatment) %>%
      summarize(
        Mean = mean(Value),
        SD = sd(Value)) %>%
      ungroup() %>%
      dplyr::select(Treatment, Trial, Mean, SD)
    data.frame()
    # Getting the plot
    plot <- ggplot(plot_df, aes(x = Trial, y = Mean, color = Treatment)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = fill_colors) +
      geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Treatment), color = NA, alpha = 0.08) + 
      scale_fill_manual(values = fill_colors) +
      theme_classic() +
      labs(y = label, x = "Trial") +
      xlim(min(plot_df$Trial), max(plot_df$Trial)) +
      scale_y_continuous(limits = c(0, NA)) +
      theme(plot.margin = margin(0.5, 0.5, 1, 0.5, "mm"),
      axis.title = element_text(size = 12, family = "Times"),
      axis.text = element_text(size = 10, family = "Times"),
      legend.position = "none")

  } else { # Plots for the intercepts and slopes for the behaviours
    if(type == "int"){
      viol_df <- df %>%
        dplyr::select(int_Control_Cold, int_CORT_Cold, int_Control_Hot, int_CORT_Hot) %>%
        pivot_longer(cols = everything(), names_to = "treatment", values_to = "values") %>%
        mutate(treatment = factor(treatment,
          levels = c("int_CORT_Hot", "int_Control_Hot", "int_CORT_Cold", "int_Control_Cold"),
          labels = c("int_CORT-Hot" = paste0("CORT-Hot (n=", n_CORT_Hot, ")"),
                    "int_Control-Hot" = paste0("Control-Hot (n=", n_Control_Hot, ")"),
                    "int_CORT-Cold" = paste0("CORT-Cold (n=", n_CORT_Cold, ")"),
                    "int_Control-Cold" = paste0("Control-Cold (n=", n_Control_Cold, ")"))))
      # Getting the real data (depending on the distribution of the response variable)
      if(res == "emergence"){
        viol_df <- viol_df %>%
          mutate(values = plogis(values))  # Convert log-odds to probabilities
      } else {
        viol_df <- viol_df %>%
          mutate(values = exp(values))  # Convert log-transformed values back to original scale
      }
    } else if(type == "slop"){
    viol_df <- df %>%
      dplyr::select(slope_Control_Cold, slope_CORT_Cold, slope_Control_Hot, slope_CORT_Hot) %>%
      pivot_longer(cols = everything(), names_to = "treatment", values_to = "values") %>%
      mutate(treatment = factor(treatment,
          levels = c("slope_CORT_Hot", "slope_Control_Hot", "slope_CORT_Cold", "slope_Control_Cold"),
          labels = c("slope_CORT_Hot" = paste0("CORT-Hot (n=", n_CORT_Hot, ")"),
                    "slope_Control_Hot" = paste0("Control-Hot (n=", n_Control_Hot, ")"),
                    "slope_CORT_Cold" = paste0("CORT-Cold (n=", n_CORT_Cold, ")"),
                    "slope_Control_Cold" = paste0("Control-Cold (n=", n_Control_Cold, ")"))))
    }

    # Getting the database for the segment plots
    segment_df <- viol_df %>%
      group_by(treatment) %>%
      summarize(
        Mean = mean(values),
        quant_min = quantile(values, 0.025, na.rm = TRUE),
        quant_max = quantile(values, 0.975, na.rm = TRUE),
        .groups = "drop") %>%
      ungroup() %>%
    data.frame()
    
    # Make the plot
    plot <- ggplot(viol_df, aes(x = treatment, y = values, fill = treatment)) +
      geom_flat_violin(alpha = 0.5) +
      scale_fill_manual(values = fill_colors) +
      geom_segment(data = segment_df, aes(y = quant_min, yend = quant_max, x = treatment, xend = treatment), size = 1.5, color = "black") +
      geom_point(data = segment_df, aes(y = Mean, x = treatment), color = "black", fill = "black", size = 3) +
      ylim(min(viol_df$values) - 0.1*mean(viol_df$values), max(viol_df$values)) +
      coord_flip() +
      guides(fill = guide_legend(reverse = TRUE)) +
      theme_classic() +
      labs(y = label, x = "Treatments") +
      theme(plot.margin = margin(0.5, 0.5, 1, 0.5, "mm"),
        axis.text.y = element_blank(),   # Remove y-axis labels
        axis.title = element_text(size = 12, family = "Times"),
        axis.text = element_text(size = 10, family = "Times"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 11, family = "Times"))
  }
  # Return the plot
  return(plot)
}