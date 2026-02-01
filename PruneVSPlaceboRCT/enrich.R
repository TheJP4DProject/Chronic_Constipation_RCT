library(clusterProfiler)
library(tidyr)
library(tidyverse)
mas_all_metadata=as.data.frame(read.csv("tableS7_KO.csv"))


conditions <- c("Placebo_04w", "Prune_04w", "Placebo_08w", "Prune_08w", "Placebo_48w", "Prune_48w")

wide_result <- data.frame(
  KO_list = character(),
  category = character(), 
  subcategory = character(),
  ID = character(),
  Description = character(),
  stringsAsFactors = FALSE
)

long_result <- data.frame()


for(cond in conditions) {
  wide_result[[paste0(cond, "_GeneRatio")]] <- character()
  wide_result[[paste0(cond, "_BgRatio")]] <- character()
  wide_result[[paste0(cond, "_RichFactor")]] <- numeric()
  wide_result[[paste0(cond, "_pvalue")]] <- numeric()
  wide_result[[paste0(cond, "_fdr")]] <- numeric()
}


for(cond in conditions) {

  pval_col <- paste0(cond, "_pval")
  fdr_col <- paste0(cond, "_fdr") 
  log2fc_col <- paste0(cond, "_log2fc")
  

  KO_up <- mas_all_metadata %>% 
    filter(!!sym(pval_col) < 0.05, !!sym(log2fc_col) > 0) %>% 
    pull(KO)
  
  KO_down <- mas_all_metadata %>% 
    filter(!!sym(pval_col) < 0.05, !!sym(log2fc_col) < 0) %>% 
    pull(KO)
  

  if(length(KO_up) > 0) {
    ko_up_en <- try({
      enrichKEGG(gene = KO_up,
                 organism = "ko", 
                 keyType = "kegg",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 1)@result %>%
        as_tibble() %>%
        mutate(metadata = cond, KO_list = "Enrich", .before = 1)
    }, silent = TRUE)
    
    if(!inherits(ko_up_en, "try-error") && nrow(ko_up_en) > 0) {
      long_result <- rbind(long_result, ko_up_en)
      
      for(j in 1:nrow(ko_up_en)) {
        pathway_id <- ko_up_en$ID[j]
        
        existing_row <- which(wide_result$ID == pathway_id & wide_result$KO_list == "Enrich")
        
        if(length(existing_row) == 0) {
          new_row <- data.frame(
            KO_list = "Enrich",
            category = ko_up_en$category[j],
            subcategory = ko_up_en$subcategory[j], 
            ID = pathway_id,
            Description = ko_up_en$Description[j],
            stringsAsFactors = FALSE
          )
          
          for(other_cond in conditions) {
            new_row[[paste0(other_cond, "_GeneRatio")]] <- NA
            new_row[[paste0(other_cond, "_BgRatio")]] <- NA
            new_row[[paste0(other_cond, "_RichFactor")]] <- NA
            new_row[[paste0(other_cond, "_pvalue")]] <- NA
            new_row[[paste0(other_cond, "_p.adjust")]] <- NA
          }
          
          wide_result <- rbind(wide_result, new_row)
          existing_row <- nrow(wide_result)
        }
        
        wide_result[existing_row, paste0(cond, "_GeneRatio")] <- ko_up_en$GeneRatio[j]
        wide_result[existing_row, paste0(cond, "_BgRatio")] <- ko_up_en$BgRatio[j]
        wide_result[existing_row, paste0(cond, "_RichFactor")] <- ko_up_en$RichFactor[j]
        wide_result[existing_row, paste0(cond, "_pvalue")] <- ko_up_en$pvalue[j]
        wide_result[existing_row, paste0(cond, "_p.adjust")] <- ko_up_en$p.adjust[j]
      }
    }
  }
  
  if(length(KO_down) > 0) {
    ko_down_en <- try({
      enrichKEGG(gene = KO_down,
                 organism = "ko",
                 keyType = "kegg", 
                 pAdjustMethod = "BH",
                 pvalueCutoff = 1)@result %>%
        as_tibble() %>%
        mutate(metadata = cond, KO_list = "Deplete", .before = 1)
    }, silent = TRUE)
    
    if(!inherits(ko_down_en, "try-error") && nrow(ko_down_en) > 0) {
      long_result <- rbind(long_result, ko_down_en)
      
      for(j in 1:nrow(ko_down_en)) {
        pathway_id <- ko_down_en$ID[j]
        

        existing_row <- which(wide_result$ID == pathway_id & wide_result$KO_list == "Deplete")
        
        if(length(existing_row) == 0) {
          new_row <- data.frame(
            KO_list = "Deplete",
            category = ko_down_en$category[j],
            subcategory = ko_down_en$subcategory[j],
            ID = pathway_id, 
            Description = ko_down_en$Description[j],
            stringsAsFactors = FALSE
          )
          
          for(other_cond in conditions) {
            new_row[[paste0(other_cond, "_GeneRatio")]] <- NA
            new_row[[paste0(other_cond, "_BgRatio")]] <- NA
            new_row[[paste0(other_cond, "_RichFactor")]] <- NA
            new_row[[paste0(other_cond, "_pvalue")]] <- NA
            new_row[[paste0(other_cond, "_p.adjust")]] <- NA
          }
          
          wide_result <- rbind(wide_result, new_row)
          existing_row <- nrow(wide_result)
        }
        
        wide_result[existing_row, paste0(cond, "_GeneRatio")] <- ko_down_en$GeneRatio[j]
        wide_result[existing_row, paste0(cond, "_BgRatio")] <- ko_down_en$BgRatio[j]
        wide_result[existing_row, paste0(cond, "_RichFactor")] <- ko_down_en$RichFactor[j]
        wide_result[existing_row, paste0(cond, "_pvalue")] <- ko_down_en$pvalue[j]
        wide_result[existing_row, paste0(cond, "_p.adjust")] <- ko_down_en$p.adjust[j]
      }
    }
  }
  
}
write.csv(wide_result, "tableS9_RCT_enrich", row.names = FALSE)






