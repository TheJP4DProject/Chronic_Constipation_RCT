setwd("/Users/gizen/Documents/prune/RCT1287/alpha_beta/")
library(openxlsx)
library(dplyr)
library(ggplot2)

meta5137=as.data.frame(read.csv("metadata.csv"))



diet_lifestyle_vars <- c("Vegetable", "Fruit", "Meat", "Processed.meat", "Seafood", 
                        "Rice", "Bread", "Noodle", "Yogurt", "Milk", 
                         "Coffee", "Alcohol", "CurrentSmoking", "Moderate.active")

adjust_vars <- c("Age", "Sex", "BMI")

existing_adjust_vars <- adjust_vars[adjust_vars %in% colnames(meta5137)]

perform_logistic_regression <- function(predictor, outcome = "Chronic.contipation", adjust_vars = NULL, data) {
  if(is.null(adjust_vars)) {
    formula_str <- paste(outcome, "~", predictor)
    analysis_type <- "Crude"
  } else {
    formula_str <- paste(outcome, "~", predictor, "+", paste(adjust_vars, collapse = " + "))
    analysis_type <- "Adjusted"
  }
  
  if(is.null(adjust_vars)) {
    reg_data <- data[, c(outcome, predictor)]
  } else {
    reg_data <- data[, c(outcome, predictor, adjust_vars)]
  }
  reg_data <- reg_data[complete.cases(reg_data), ]
  
  if(nrow(reg_data) == 0) {
    return(list(or = NA, pvalue = NA, ci_lower = NA, ci_upper = NA, n = 0, analysis_type = analysis_type))
  }
  
  model <- glm(as.formula(formula_str), data = reg_data, family = binomial)
  coef_summary <- summary(model)$coefficients
  
  ci <- confint(model)
  
  coef <- coef_summary[predictor, "Estimate"]
  or <- exp(coef)
  ci_lower <- exp(ci[predictor, 1])
  ci_upper <- exp(ci[predictor, 2])
  pvalue <- coef_summary[predictor, "Pr(>|z|)"]
  
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

results_list <- list()

for(predictor in existing_vars) {
  crude_result <- perform_logistic_regression(predictor, "Chronic.contipation", adjust_vars = NULL, meta5137)
  results_list[[paste0(predictor, "_crude")]] <- crude_result
  adj_result <- perform_logistic_regression(predictor, "Chronic.contipation", existing_adjust_vars, meta5137)
  results_list[[paste0(predictor, "_adj")]] <- adj_result
}

# 准备森林图数据
create_forest_data <- function(results_list, predictors) {
  forest_data <- data.frame()
  
  for(predictor in predictors) {
    crude_key <- paste0(predictor, "_crude")
    adj_key <- paste0(predictor, "_adj")
    
    # 跳过不存在的结果
    if(!crude_key %in% names(results_list) || !adj_key %in% names(results_list)) next
    
    # Crude结果
    crude_row <- data.frame(
      Predictor = predictor,
      Analysis = "Crude",
      OR = results_list[[crude_key]]$or,
      CI_lower = results_list[[crude_key]]$ci_lower,
      CI_upper = results_list[[crude_key]]$ci_upper,
      P_value = results_list[[crude_key]]$pvalue,
      N = results_list[[crude_key]]$n
    )
    
    # Adjusted结果
    adj_row <- data.frame(
      Predictor = predictor,
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

forest_data <- create_forest_data(results_list, existing_vars)

write.csv(forest_data, "tableS1_PART1.csv", row.names = FALSE)


diet_lifestyle_vars <- c(
  "Vegetable", "Fruit", "Meat", "Processed.meat", "Seafood",
  "Rice", "Bread", "Noodle", "Yogurt", "Milk",
  "Coffee", "Alcohol", "CurrentSmoking", "Moderate.active",
  "Age", "Sex", "BMI"
)


df <- meta5137 %>%
  select(Chronic.contipation, all_of(diet_lifestyle_vars))


is_binary <- function(x) {
  ux <- unique(na.omit(x))
  all(ux %in% c(0, 1))
}


summarize_var <- function(data, var) {
  x <- data[[var]]
  
  if (is_binary(x)) {
    prop <- mean(x == 1, na.rm = TRUE) * 100
    return(sprintf("%.1f%%", prop))
  } else {
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    return(sprintf("%.2f (%.2f)", m, s))
  }
}

result <- data.frame(
  Variable = diet_lifestyle_vars,
  health = sapply(diet_lifestyle_vars, function(v)
    summarize_var(df %>% filter(Chronic.contipation == 0), v)),
  Chronic.contipation = sapply(diet_lifestyle_vars, function(v)
    summarize_var(df %>% filter(Chronic.contipation == 1), v)),
  stringsAsFactors = FALSE
)
result
write.csv(result,"tableS1_PART2.csv")
