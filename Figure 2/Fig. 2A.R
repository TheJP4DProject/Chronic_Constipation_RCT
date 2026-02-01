library(tidyverse)
library(openxlsx)  
library(vegan)    



meta5137=as.data.frame(read.csv("metadata.csv"))
fdata=as.data.frame(read.csv("motusp.csv"))


set.seed(123)
data=fdata[rownames(meta5137),]

identical(rownames(data),rownames(meta5137))


data_minprev <- data[,(colSums(data>0)/nrow(data))>0.025]

meta5137$shannon_minprev <- diversity(data_minprev, index = "shannon")
meta5137$richness_minprev <- scale(log10(rowSums(data_minprev > 0)))

ml=read.table("motusp_ml.tsv",sep="\t",header = 1,row.names = 1)
modify_rownames <- function(df) {
  old_names <- rownames(df)
  new_names <- sapply(old_names, function(name) {
    num_name <- suppressWarnings(as.numeric(name))
    if (!is.na(num_name)) {
      return(as.character(num_name))
    } else {
      return(name)
    }
  })
  rownames(df) <- new_names
  return(df)
}

ml=modify_rownames(ml)
meta5137$ml=ml[meta5137$METAF,]


meta5137$ml_scaled <- scale(log10(meta5137$ml))[,1]
adjust_vars <- colnames(meta5137)[4:12]  # Age, Sex, BMI, Vegetable, Fruit, Yogurt, Bread, Alcohol, Moderate.active

perform_regression <- function(outcome_var, predictor = "Chronic.contipation", adjust_vars = NULL, data) {
  if(is.null(adjust_vars)) {
    formula_str <- paste(outcome_var, "~", predictor)
    analysis_type <- "Crude"
  } else {
    formula_str <- paste(outcome_var, "~", predictor, "+", paste(adjust_vars, collapse = " + "))
    analysis_type <- "Adjusted"
  }
  
  if(is.null(adjust_vars)) {
    reg_data <- data[, c(outcome_var, predictor)]
  } else {
    reg_data <- data[, c(outcome_var, predictor, adjust_vars)]
  }
  reg_data <- reg_data[complete.cases(reg_data), ]
  
  if(nrow(reg_data) == 0) {
    return(list(coef = NA, pvalue = NA, ci_lower = NA, ci_upper = NA, n = 0, analysis_type = analysis_type))
  }
  
  model <- lm(as.formula(formula_str), data = reg_data)
  coef_summary <- summary(model)$coefficients
  
  ci <- confint(model)
  coef <- coef_summary[predictor, "Estimate"]
  or <- exp(coef)
  ci_lower <- exp(ci[predictor, 1])
  ci_upper <- exp(ci[predictor, 2])
  pvalue <- coef_summary[predictor, "Pr(>|t|)"]
  
  return(list(
    or = or,
    coef = coef,
    pvalue = pvalue,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    n = nrow(reg_data),
    analysis_type = analysis_type
  ))
}

outcomes <- c("shannon_minprev", "richness_minprev", "ml_scaled")

results_list <- list()

for(outcome in outcomes) {
  crude_result <- perform_regression(outcome, "Chronic.contipation", adjust_vars = NULL, meta5137)
  results_list[[paste0(outcome, "_crude")]] <- crude_result
  
  adj_result <- perform_regression(outcome, "Chronic.contipation", adjust_vars, meta5137)
  results_list[[paste0(outcome, "_adj")]] <- adj_result
}

create_forest_data <- function(results_list, outcomes) {
  forest_data <- data.frame()
  
  for(outcome in outcomes) {
    crude_key <- paste0(outcome, "_crude")
    adj_key <- paste0(outcome, "_adj")
    
    if(outcome == "shannon_minprev") {
      category <- "Shannon Diversity"
    } else if(outcome == "richness_minprev") {
      category <- "Richness"
    } else if(outcome == "ml_scaled") {
      category <- "Microbial Load"
    }
    
    # Crude结果
    crude_row <- data.frame(
      Category = category,
      Analysis = "Crude",
      OR = results_list[[crude_key]]$or,
      CI_lower = results_list[[crude_key]]$ci_lower,
      CI_upper = results_list[[crude_key]]$ci_upper,
      P_value = results_list[[crude_key]]$pvalue,
      N = results_list[[crude_key]]$n
    )
    
    # Adjusted结果
    adj_row <- data.frame(
      Category = category,
      Analysis = "Adjusted",
      OR = results_list[[adj_key]]$or,
      CI_lower = results_list[[adj_key]]$ci_lower,
      CI_upper = results_list[[adj_key]]$ci_upper,
      P_value = results_list[[adj_key]]$pvalue,
      N = results_list[[adj_key]]$n
    )
    
    forest_data <- rbind(forest_data, crude_row, adj_row)
  }
  
  return(forest_data)
}

forest_data <- create_forest_data(results_list, outcomes)

colors <- c("Crude" = "#2E86AB", "Adjusted" = "#A23B72")

forest_data$Category <- factor(forest_data$Category, 
                               levels = c("Shannon Diversity", "Richness", "Microbial Load"))
forest_data$Analysis <- factor(forest_data$Analysis, 
                               levels = c( "Adjusted","Crude"))

forest_data$Significance <- ifelse(forest_data$P_value < 0.001, "***",
                                   ifelse(forest_data$P_value < 0.01, "**",
                                          ifelse(forest_data$P_value < 0.05, "*", "ns")))

forest_data$Label <- paste0("OR: ", sprintf("%.3f", forest_data$OR), 
                            "\n95% CI: (", sprintf("%.3f", forest_data$CI_lower), 
                            ", ", sprintf("%.3f", forest_data$CI_upper), ")",
                            "\nP: ", sprintf("%.3f", forest_data$P_value),
                            "\nN: ", forest_data$N)

p <- ggplot(forest_data, aes(x = Category, y = OR, color = Analysis)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                position = position_dodge(width = 0.5), width = 0.2) +
  scale_color_manual(values = colors) +
  scale_y_continuous(breaks = function(x) {
    range_vals <- range(x, na.rm = TRUE)
    breaks <- pretty(range_vals, n = 6)
    if(!1.0 %in% breaks) {
      breaks <- sort(c(breaks, 1.0))
    }
    return(breaks)
  }, labels = function(x) {
    sapply(x, function(val) {
      if(is.na(val)) {
        return(as.character(val))
      } else if(abs(val - 1.0) < 0.001) {
        return("1.0")
      } else {
        return(as.character(val))
      }
    })
  }) +
  labs(title = "Chronic Constipation",
       subtitle = "MinPrev≥0.025 | Adjusted: Age, Sex, BMI, Vegetable, Fruit, Yogurt, Bread, Alcohol, Moderate activity",
       x = "Alpha Diversity Index",
       y = "Odds Ratio (95% CI)",
       color = "Analysis Type") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        panel.grid.minor = element_blank()) +
  coord_flip()
print(p)

ggsave("Fig2A.pdf", p, width = 4.5, height = 3.5, dpi = 300)
