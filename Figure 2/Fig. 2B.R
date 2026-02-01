library(vegan)
library(ggplot2)
library(gridExtra)

meta5137=as.data.frame(read.csv("metadata.csv"))
fdata=as.data.frame(read.csv("motusp.csv"))

set.seed(123)
data=fdata[rownames(meta5137),]

data_minprev <- data[,(colSums(data>0)/nrow(data))>0.025]
adjust_vars <- colnames(meta5137)[4:12]

meta_complete <- meta5137[complete.cases(meta5137[, c("Chronic.contipation", adjust_vars)]), ]
data_minprev_complete <- data_minprev[rownames(meta_complete), ] / rowSums(data_minprev[rownames(meta_complete), ])

bc_dist_minprev <- vegdist(data_minprev_complete, method = "bray")

formula_vars <- c("Chronic.contipation", adjust_vars)
formula_str <- paste("dist ~", paste(formula_vars, collapse = " + "))

adonis_minprev <- adonis2(bc_dist_minprev ~ ., data = meta_complete[, formula_vars], permutations = 999, by = "terms")

cc_minprev <- adonis_minprev[1, ]

mds_minprev <- metaMDS(bc_dist_minprev, k = 2, trymax = 20)

create_plot <- function(ord_data, metadata, title, r2, pval, method = "MDS") {
  if(method == "MDS") {
    plot_data <- data.frame(
      x = ord_data$points[, 1],
      y = ord_data$points[, 2],
      group = factor(metadata$Chronic.contipation, levels = c(0, 1), labels = c("Control", "Constipation"))
    )
    x_label <- "MDS1"
    y_label <- "MDS2"
    stats_text <- paste0("R² = ", sprintf("%.4f", r2), "\nP = ", sprintf("%.3f", pval))
  } else {
    explained <- ord_data$values$Relative_eig[1:2] * 100
    plot_data <- data.frame(
      x = ord_data$vectors[, 1],
      y = ord_data$vectors[, 2],
      group = factor(metadata$Chronic.contipation, levels = c(0, 1), labels = c("Control", "Constipation"))
    )
    x_label <- paste0("PCo1 (", round(explained[1], 1), "%)")
    y_label <- paste0("PCo2 (", round(explained[2], 1), "%)")
    stats_text <- paste0("R² = ", sprintf("%.4f", r2), "\nP = ", sprintf("%.3f", pval))
  }
  
  ggplot(plot_data, aes(x = x, y = y, color = group)) +
    geom_point(size = 2, alpha = 0.7) +
    stat_ellipse(type = "norm", level = 0.95, size = 1) +
    scale_color_manual(values = c("Control" = "#2E86AB", "Constipation" = "#A23B72")) +
    labs(title = title, x = x_label, y = y_label) +
    theme_minimal() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold")) +
    annotate("text", x = -Inf, y = Inf, label = stats_text, hjust = -0.1, vjust = 1.2, 
             size = 3.5, fontface = "bold")
}

p1 <- create_plot(mds_minprev, meta_complete, "MDS - MinPrev≥0.025", cc_minprev$R2, cc_minprev$`Pr(>F)`, "MDS")

pdf("Fig2B.pdf", width = 6, height = 5)
grid.arrange(
  p1,
  ncol = 1,
  top = "Beta Diversity Analysis - Chronic Constipation"
)
dev.off()
