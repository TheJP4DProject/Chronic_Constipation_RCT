library(ggplot2)
library(dplyr)



maa <- read.csv("tableS2_motusp.csv") 
wilcox <-  read.csv("tableS3_motusp.csv") 

maa$GTDB=gsub("\\[ID:(\\d+)\\]", "\\1", maa$GTDB)
wilcox$GTDB=gsub("\\[ID:(\\d+)\\]", "\\1", wilcox$GTDB)
maa=maa[maa$qval<0.1,]
wilcox$Placebo_04w_pval=ifelse(is.na(wilcox$Placebo_04w_pval),1,wilcox$Placebo_04w_pval)
wilcox$Placebo_08w_pval=ifelse(is.na(wilcox$Placebo_08w_pval),1,wilcox$Placebo_08w_pval)
wilcox$Prune_04w_pval=ifelse(is.na(wilcox$Prune_04w_pval),1,wilcox$Prune_04w_pval)
wilcox$Prune_08w_pval=ifelse(is.na(wilcox$Prune_08w_pval),1,wilcox$Prune_08w_pval)
wilcox=wilcox[(wilcox$Placebo_04w_pval<0.05)|(wilcox$Placebo_08w_pval<0.05)|(wilcox$Prune_04w_pval<0.05)|(wilcox$Prune_08w_pval<0.05),]
maa=maa[maa$GTDB %in% wilcox$GTDB,]
wilcox=wilcox[wilcox$GTDB %in% maa$GTDB,]


rownames(maa)=maa$GTDB
rownames(wilcox)=wilcox$GTDB
maa=maa[rownames(wilcox),]
lg2fc=wilcox[,c(7,13,4,10)]
p=wilcox[,c(8,14,5,11)]
lg2fc[is.na(lg2fc)]=0
lg2fc[p>0.05]=0

lg2fc[lg2fc>0]=1
lg2fc[lg2fc<0]=-1
lg2fc <- apply(lg2fc, 2, as.numeric)
col_rho_1 = colorRamp2(
  c(-1, 0,1), 
  c("lightblue", "white", "salmon")
)

hm2 <- Heatmap(t(as.matrix(lg2fc)), 
               name = "lg2fc",
               cluster_rows = FALSE,
               cluster_columns = F,
               row_names_side = "left",show_column_dend = F,
               height = nrow(t(lg2fc)) * unit(2, "mm"),
               width = ncol(t(lg2fc)) * unit(2, "mm"),
               rect_gp = gpar(col = "gray", lwd =  0.5),
               column_names_gp = gpar(fontsize = 5, fontface="bold"),
               row_names_gp = gpar(fontsize = 5, fontface="bold"),
               col = col_rho_1,
               column_title_gp = gpar(fontsize = 15),
               heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface="bold"), labels_gp = gpar(fontsize = 10),
                                           legend_direction = "horizontal"))

col_rho_2 = colorRamp2(
  c(-0.5 * max(abs(maa$coef)), 0, 0.5 * max(abs(maa$coef))), 
  c("blue", "white", "red")
)

hm1 <- Heatmap(t(as.matrix(maa[,"coef",drop=F])), 
               name = "*FDR<0.05\n**FDR<0.01\nMAA coef",
               cluster_rows = FALSE,
               cluster_columns = T,
               row_names_side = "left",show_column_dend = F,
               height = nrow(t(maa[,"coef",drop=F])) * unit(2, "mm"),
               width = ncol(t(maa[,"coef",drop=F])) * unit(2, "mm"),
               rect_gp = gpar(col = "gray", lwd = 0.5),
               column_names_gp = gpar(fontsize = 5, fontface="bold"),
               row_names_gp = gpar(fontsize = 5, fontface="bold"),
               col = col_rho_2,
               column_title_gp = gpar(fontsize = 15),
               heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface="bold"), labels_gp = gpar(fontsize = 10),
                                           legend_direction = "horizontal", at = c(-0.1, -0.05, 0, 0.05, 0.1)),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(t(maa[,"qval",drop=F])[i, j] < 0.05) {
                   if(t(maa[,"qval",drop=F])[i, j] < 0.01){grid.text("**", x, y, vjust = 0.8,gp = gpar(fontsize = 5))}
                   else{grid.text("*", x, y, vjust = 0.8,gp = gpar(fontsize = 5))}
                   }})
set.seed(123)
cairo_pdf(file = paste0( "Fig3D.pdf"),  family="Helvetica",width = 10)
ht5 = draw(hm2 %v% hm1, padding = unit(c(80, 50, 2, 2), "mm"), merge_legend=TRUE,main="*FDR<0.05\n**FDR<0.01\nMAA coef"
)
dev.off()

