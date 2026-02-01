
library(openxlsx)
library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)

meta5137 <- as.data.frame(read.csv("metadata.csv"))

fdata <- as.data.frame(read.csv("motusp.csv"))
available_samples <- rownames(fdata)

prepare_groups <- function(prune_group){
  df0 <- meta5137 %>% filter(Prune==prune_group & !is.na(`0WMetaf`)) %>%
    select(ID, Sample=`0WMetaf`) %>% mutate(Time="0W")
  df4 <- meta5137 %>% filter(Prune==prune_group) %>%
    select(ID, Sample=MetaF) %>% mutate(Time="4W")
  df8 <- meta5137 %>% filter(Prune==prune_group & !is.na(`8WMETAF`)) %>%
    select(ID, Sample=`8WMETAF`) %>% mutate(Time="8W")
  
  groups_df <- bind_rows(df0, df4, df8) %>% filter(Sample %in% available_samples)
  return(groups_df)
}
perform_adonis_strata <- function(fdata_matrix, groups_df, comparison_name){
  # Keep subjects with samples at both time points
  subjects_keep <- groups_df %>%
    group_by(ID) %>%
    summarise(n_time = n_distinct(Time)) %>%
    filter(n_time == 2) %>% pull(ID)
  
  groups_df <- groups_df %>% filter(ID %in% subjects_keep)
  fdata_matrix <- fdata_matrix[groups_df$Sample, ]
  
  if(length(unique(groups_df$Time)) < 2){
    cat(comparison_name, ": Only one timepoint left, skipping\n")
    return(NULL)
  }
  
  result <- tryCatch({
    adonis2(fdata_matrix ~ Time, data = groups_df, method="bray",
            permutations=999, strata=groups_df$ID)
  }, error=function(e){
    cat("ERROR in", comparison_name, "\n")
    return(NULL)
  })
  
  return(result)
}

format_adonis_result <- function(adonis_result){
  if(is.null(adonis_result)) return(list(r2=NA, p_text="N/A"))
  r2 <- round(adonis_result$R2[1],3)
  p <- adonis_result$`Pr(>F)`[1]
  p_text <- ifelse(p<0.001, "p < 0.001", paste0("p = ", round(p,3)))
  list(r2=r2, p_text=p_text)
}

plot_ordination <- function(fdata_matrix, groups_df, title_prefix, colors){
  # Distance
  bc_dist <- vegdist(fdata_matrix, method="bray")
  
  # NMDS
  set.seed(123)
  nmds <- metaMDS(fdata_matrix, distance="bray", k=2, trymax=100)
  
  # PERMANOVA comparisons
  comparisons <- list(c("0W","4W"), c("0W","8W"))
  labels <- c()
  for(comp in comparisons){
    idx <- groups_df$Time %in% comp
    ad <- perform_adonis_strata(fdata_matrix[idx,], groups_df[idx,], paste0(title_prefix, ": ", comp[1]," vs ",comp[2]))
    res <- format_adonis_result(ad)
    if(!is.na(res$r2)){
      labels <- c(labels, paste0(comp[1]," vs ",comp[2],": R²=",res$r2,", ",res$p_text))
    }
  }
  label_text <- paste(labels, collapse="\n")
  
  
  # NMDS plot
  plot_nmds <- data.frame(NMDS1=nmds$points[,1], NMDS2=nmds$points[,2], Group=groups_df$Time)
  p_nmds <- ggplot(plot_nmds, aes(NMDS1, NMDS2, color=Group)) +
    stat_ellipse(level=0.95, size=0.5) +
    geom_point(size=1, alpha=0.7) +
    scale_color_manual(values=colors) +
    labs(x="NMDS1", y="NMDS2",
         title=paste0("NMDS - ", title_prefix," (Stress=", round(nmds$stress,3),")")) +
    theme_bw(base_family = "Helvetica", base_size = 5)
  if(length(labels)>0){
    p_nmds <- p_nmds + annotate("text", x=-Inf, y=Inf, label=label_text,
                                hjust=-0.05, vjust=1.1,
                                size = 5 / 2.845,     
                                family = "Helvetica",
                                fontface = "bold",
                                lineheight = 0.9)
  }
  
  return(list(pcoa=p_pcoa, nmds=p_nmds))
}

groups_P1 <- prepare_groups(1)
fdata_P1 <- fdata[groups_P1$Sample, ]
fdata_P1 <- fdata_P1[, colSums(fdata_P1)>0]
colors_P1 <- c("0W"="#E6E6FA","4W"="#9370DB","8W"="#8A2BE2")
plots_P1 <- plot_ordination(fdata_P1, groups_P1, "Prune", colors_P1)

groups_P0 <- prepare_groups(0)
fdata_P0 <- fdata[groups_P0$Sample, ]
fdata_P0 <- fdata_P0[, colSums(fdata_P0)>0]
colors_P0 <- c("0W"="#87CEEB","4W"="#4682B4","8W"="#1E90FF")

plots_P0 <- plot_ordination(fdata_P0, groups_P0, "Placebo", colors_P0)

pdf("PCoA_NMDS_Prune.pdf", width=4.5, height=1.8)
grid.arrange( plots_P1$nmds,plots_P0$nmds, ncol=2)
dev.off()

