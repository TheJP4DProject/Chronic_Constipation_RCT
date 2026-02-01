library(tidyverse)   
library(vegan)       
library(gridExtra)   

meta5137=as.data.frame(read.csv("metadata.csv"))
meta5137=meta5137[meta5137$Chronic.contipation==1,]

fdata=as.data.frame(read.csv("motusp.csv"))

data=fdata
data=data[,(colSums(data>0)/nrow(data))>0]



prune0w=as.data.frame(read.csv("RCT_0wid.csv"))
prune4w=as.data.frame(read.csv("RCT_4wid.csv"))
prune8w=as.data.frame(read.csv("RCT_8wid.csv"))
rownames(prune0w)=prune0w$ID
rownames(prune4w)=prune4w$ID
rownames(prune8w)=prune8w$ID

# 定义函数计算alpha多样性指标
calculate_alpha_diversity <- function(data_matrix) {
  # Shannon多样性
  shannon <- diversity(data_matrix, index = "shannon")
  # Richness (log10转换)
  richness <- rowSums(data_matrix > 0)
  
  return(list(shannon = shannon, richness = richness))
}

prune0w1key=prune0w[prune0w$Prune==1,]$Metaf
prune0w0key=prune0w[prune0w$Prune==0,]$Metaf
prune4w1key=prune4w[prune4w$Prune==1,]$MetaF
prune4w0key=prune4w[prune4w$Prune==0,]$MetaF
prune8w1key=prune8w[prune8w$Prune==1,]$METAF
prune8w0key=prune8w[prune8w$Prune==0,]$METAF

prune0w1key=prune0w1key[prune0w1key %in% rownames(fdata)]
prune0w0key=prune0w0key[prune0w0key %in% rownames(fdata)]
prune4w1key=prune4w1key[prune4w1key %in% rownames(fdata)]
prune4w0key=prune4w0key[prune4w0key %in% rownames(fdata)]
prune8w1key=prune8w1key[prune8w1key %in% rownames(fdata)]
prune8w0key=prune8w0key[prune8w0key %in% rownames(fdata)]

alpha_0w_prune1 <- calculate_alpha_diversity(data[prune0w1key,])
alpha_0w_prune0 <- calculate_alpha_diversity(data[prune0w0key,])
alpha_4w_prune1 <- calculate_alpha_diversity(data[prune4w1key,])
alpha_4w_prune0 <- calculate_alpha_diversity(data[prune4w0key,])
alpha_8w_prune1 <- calculate_alpha_diversity(data[prune8w1key,])
alpha_8w_prune0 <- calculate_alpha_diversity(data[prune8w0key,])

# 创建数据框用于绘图 - 使用所有数据
create_plot_data <- function(data_0w_p0, data_0w_p1, data_4w_p0, data_4w_p1, data_8w_p0, data_8w_p1, metric_name) {
  plot_data <- data.frame(
    Value = c(data_0w_p0, data_0w_p1, data_4w_p0, data_4w_p1, data_8w_p0, data_8w_p1),
    Time = c(rep("0w", length(data_0w_p0) + length(data_0w_p1)),
             rep("4w", length(data_4w_p0) + length(data_4w_p1)),
             rep("8w", length(data_8w_p0) + length(data_8w_p1))),
    Group = c(rep("Placebo", length(data_0w_p0)), rep("Prune", length(data_0w_p1)),
              rep("Placebo", length(data_4w_p0)), rep("Prune", length(data_4w_p1)),
              rep("Placebo", length(data_8w_p0)), rep("Prune", length(data_8w_p1))),
    Metric = metric_name
  )
  return(plot_data)
}

shannon_data <- create_plot_data(alpha_0w_prune0$shannon, alpha_0w_prune1$shannon,
                                 alpha_4w_prune0$shannon, alpha_4w_prune1$shannon,
                                 alpha_8w_prune0$shannon, alpha_8w_prune1$shannon, "Shannon")




perform_stats <- function(metric_name) {
  

  prune0w1keyid <- prune0w[(prune0w$Prune==1)&(prune0w$Metaf %in% prune0w1key),]$ID
  prune0w0keyid <- prune0w[(prune0w$Prune==0)&(prune0w$Metaf %in% prune0w0key),]$ID
  prune4w1keyid <- prune4w[(prune4w$Prune==1)&(prune4w$MetaF %in% prune4w1key),]$ID
  prune4w0keyid <- prune4w[(prune4w$Prune==0)&(prune4w$MetaF %in% prune4w0key),]$ID
  prune8w1keyid <- prune8w[(prune8w$Prune==1)&(prune8w$METAF %in% prune8w1key),]$ID
  prune8w0keyid <- prune8w[(prune8w$Prune==0)&(prune8w$METAF %in% prune8w0key),]$ID
  
  check04w1 <- intersect(prune0w1keyid, prune4w1keyid)
  check04w0 <- intersect(prune0w0keyid, prune4w0keyid)
  check08w1 <- intersect(prune0w1keyid, prune8w1keyid)
  check08w0 <- intersect(prune0w0keyid, prune8w0keyid)
  
  results <- list()

  if(metric_name == "Shannon") {
    data_0w_p1 <- alpha_0w_prune1$shannon
    data_0w_p0 <- alpha_0w_prune0$shannon
    data_4w_p1 <- alpha_4w_prune1$shannon
    data_4w_p0 <- alpha_4w_prune0$shannon
    data_8w_p1 <- alpha_8w_prune1$shannon
    data_8w_p0 <- alpha_8w_prune0$shannon
  } 
  
  wilcox_pair <- function(check, data_tw1, data_tw2, df_tw1, df_tw2, results_name, col_tw1, col_tw2) {
    if(length(check) > 0) {
      metaf_tw1_ordered <- df_tw1[df_tw1$ID %in% check, ]
      metaf_tw1_ordered <- metaf_tw1_ordered[order(metaf_tw1_ordered$ID), ][[col_tw1]]

      metaf_tw2_ordered <- df_tw2[df_tw2$ID %in% check, ]
      metaf_tw2_ordered <- metaf_tw2_ordered[order(metaf_tw2_ordered$ID), ][[col_tw2]]
      
      paired_tw1 <- data_tw1[metaf_tw1_ordered]
      paired_tw2 <- data_tw2[metaf_tw2_ordered]
      valid <- !is.na(paired_tw1) & !is.na(paired_tw2)
      paired_tw1 <- paired_tw1[valid]
      paired_tw2 <- paired_tw2[valid]
      if(length(paired_tw1) > 0) {
        results[[results_name]] <<- wilcox.test(paired_tw1, paired_tw2, paired = TRUE)$p.value
      }
    }
  }
  
  # 0w vs 4w
  wilcox_pair(check04w0, data_0w_p0, data_4w_p0, prune0w, prune4w, "placebo_04w", "Metaf", "MetaF")
  wilcox_pair(check04w1, data_0w_p1, data_4w_p1, prune0w, prune4w, "prune_04w", "Metaf", "MetaF")
  
  # 0w vs 8w
  wilcox_pair(check08w0, data_0w_p0, data_8w_p0, prune0w, prune8w, "placebo_08w", "Metaf", "METAF")
  wilcox_pair(check08w1, data_0w_p1, data_8w_p1, prune0w, prune8w, "prune_08w", "Metaf", "METAF")
  
  
  # ------------------------------
  # 4. 趋势分析（0/4/8w）使用 gls
  # ------------------------------
  build_long <- function(data_0w, data_4w, data_8w, df_0w, df_4w, df_8w) {
    sample_0w <- names(data_0w)
    sample_4w <- names(data_4w)
    sample_8w <- names(data_8w)
    
    df_long <- data.frame(
      sample = c(sample_0w, sample_4w, sample_8w),
      value  = c(data_0w, data_4w, data_8w),
      time_num = c(rep(0,length(sample_0w)),
                   rep(4,length(sample_4w)),
                   rep(8,length(sample_8w)))
    )
    
    df_long$subject <- factor(make.names(c(
      df_0w$ID[match(sample_0w, df_0w$Metaf)],
      df_4w$ID[match(sample_4w, df_4w$MetaF)],
      df_8w$ID[match(sample_8w, df_8w$METAF)]
    )))
    
    df_long <- df_long[!is.na(df_long$value), ]
    df_long$value_scaled <- as.numeric(scale(df_long$value))
    
    tab <- table(df_long$subject)
    df_long <- df_long[df_long$subject %in% names(tab[tab>=3]), ]
    
    return(df_long)
  }
  
  return(results)
}

create_boxplot <- function(data_subset, title, stats_results) {
  data_subset$Group <- factor(data_subset$Group, levels = c("Prune", "Placebo"))

  data_subset$Time_Group <- interaction(data_subset$Time, data_subset$Group, sep = "_")

  colors <- c(
    "0w_Prune" = "#E6E6FA",
    "4w_Prune" = "#9370DB", 
    "8w_Prune" = "#8A2BE2",  
    "0w_Placebo" = "#87CEEB", 
    "4w_Placebo" = "#4682B4", 
    "8w_Placebo" = "#1E90FF" 
  )
  
  p <- ggplot(data_subset, aes(x = Time, y = Value)) +
    geom_boxplot(aes(fill = Time_Group), position = position_dodge(width = 0.8), 
                 outlier.shape = NA, alpha = 0.8) +
    geom_point(aes(color = Group), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), 
               alpha = 0.6, size = 1) +
    scale_fill_manual(values = colors, guide = "none") +
    scale_color_manual(values = c("Prune" = "#9370DB", "Placebo" = "#4682B4")) + 
    facet_wrap(~ Group, ncol = 2)+
    labs(title = title, x = "Time Point") +
    theme_minimal(base_size = 10) + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold", family = "Helvetica"),  # 5pt Helvetica
      legend.position = "right",
      legend.text = element_text(size = 10, family = "Helvetica"),  # 5pt Helvetica
      legend.title = element_text(size = 10, family = "Helvetica"),  # 5pt Helvetica
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 10, face = "bold", family = "Helvetica"),  # 5pt Helvetica
      axis.text = element_text(size = 10, family = "Helvetica"),  # 5pt Helvetica
      axis.title = element_text(size = 10, family = "Helvetica"),  # 5pt Helvetica
      axis.text.x = element_text(size = 10, family = "Helvetica"),  # 5pt Helvetica
      axis.text.y = element_text(size = 10, family = "Helvetica")  # 5pt Helvetica
    )
  
  prune_label <- ""
  if(!is.null(stats_results$prune_04w)) {
    prune_label <- paste0("0w vs 4w: ", get_significance_label(stats_results$prune_04w))
  }
  if(!is.null(stats_results$prune_08w)) {
    if(prune_label != "") prune_label <- paste0(prune_label, "\n")
    prune_label <- paste0(prune_label, "0w vs 8w: ", get_significance_label(stats_results$prune_08w))
  }

  
  placebo_label <- ""
  if(!is.null(stats_results$placebo_04w)) {
    placebo_label <- paste0("0w vs 4w: ", get_significance_label(stats_results$placebo_04w))
  }
  if(!is.null(stats_results$placebo_08w)) {
    if(placebo_label != "") placebo_label <- paste0(placebo_label, "\n")
    placebo_label <- paste0(placebo_label, "0w vs 8w: ", get_significance_label(stats_results$placebo_08w))
  }


  if(placebo_label != "" || prune_label != "") {
    annotation_data <- data.frame(
      x = -Inf,
      y = Inf,
      label = c(prune_label, placebo_label),
      Group = factor(c("Prune", "Placebo"), levels = c("Prune", "Placebo")), 
      stringsAsFactors = FALSE
    )

    annotation_data <- annotation_data[annotation_data$label != "", ]
    
    if(nrow(annotation_data) > 0) {
      p <- p + geom_text(data = annotation_data, 
                         aes(x = x, y = y, label = label),
                         hjust = -0.1, vjust = 1.1, size = 5/2.835, 
                         family = "Helvetica", 
                         inherit.aes = FALSE)
    }
  }
  
  return(p)
}

shannon_stats <- perform_stats("Shannon")

get_significance_label <- function(p_value) {
  if(is.na(p_value) || is.null(p_value)) return("")
  p_text <- sprintf("p=%.3f", p_value)
  if(p_value < 0.001) return(paste0(p_text, "***"))
  if(p_value < 0.01) return(paste0(p_text, "**"))
  if(p_value < 0.05) return(paste0(p_text, "*"))
  return(p_text)
}


p1 <- create_boxplot(shannon_data, "Shannon Diversity", shannon_stats)


combined_plot <- grid.arrange(p1, ncol = 1)
ggsave("Fig3B.pdf", combined_plot, width = 5, height = 4, dpi = 300)

