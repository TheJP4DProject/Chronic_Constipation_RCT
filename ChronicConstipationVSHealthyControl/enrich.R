library(clusterProfiler)
library(tidyr)
library(tidyverse)
mas_all_metadata=as.data.frame(read.csv("tableS4_KO.csv"))

resultp=NULL

for (i in unique(mas_all_metadata$metadata)){
  temp=mas_all_metadata
  KO_up <- temp %>% filter(pval<0.05,coef>0, metadata==i) %>% pull(feature)
  KO_down <- temp %>% filter(pval<0.05,coef<0, metadata==i) %>% pull(feature)
  
  ko_up_en <- try(
    enrichKEGG(gene = KO_up,
               organism = "ko",
               keyType = "kegg",
               pAdjustMethod = "BH",
               pvalueCutoff = 1) %>%
      .@result %>% as_tibble() %>%
      mutate(metadata=i,KO_list="Enrich", .before=1)
  )
  ko_up_en$N=temp$N[1]
  ko_down_en <- try(
    enrichKEGG(gene = KO_down,
               organism = "ko",
               keyType = "kegg",
               pAdjustMethod = "BH",
               pvalueCutoff = 1) %>%
      .@result %>% as_tibble() %>%
      mutate(metadata=i,KO_list="Deplete", .before=1)
  )
  ko_down_en$N=temp$N[1]
  resultp=rbind(resultp,ko_up_en,ko_down_en)
  
}
resultp=resultp[,-17]
write.csv(resultp,"tableS6_enrich.csv")
