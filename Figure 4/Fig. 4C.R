library(dplyr)
library(ComplexHeatmap)

path <- read.csv("selected_data_from_KOandMO.csv") 

rownames(path)=paste0(path$description,", ",path$ID)
maa=path[,"coef",drop=F]
lg2fc=path[,c(10,16,7,13)]
col_rho_1 = colorRamp2(
  c(-1, 0,1), 
  c("lightblue", "white", "salmon")
)
p=path[,c(11,17,8,14)]
hm2 <- Heatmap((as.matrix(lg2fc)), 
               name = "*pval<0.05\nlg2fc",
               cluster_rows = FALSE,
               cluster_columns = F,
               row_names_side = "left",show_column_dend = F,
               height = nrow((lg2fc)) * unit(2, "mm"),
               width = ncol((lg2fc)) * unit(2, "mm"),
               rect_gp = gpar(col = "gray", lwd =  0.5),
               column_names_gp = gpar(fontsize = 5, fontface="bold"),
               row_names_gp = gpar(fontsize = 5, fontface="bold"),
               col = col_rho_1,
               column_title_gp = gpar(fontsize = 15),
               heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface="bold"), labels_gp = gpar(fontsize = 10),
                                           legend_direction = "horizontal"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if((p)[i, j] < 0.05) {grid.text("*", x, y, vjust = 0.8,gp = gpar(fontsize = 5))}})

col_rho_2 = colorRamp2(
  c(-0.5 * max(abs(maa$coef)), 0, 0.5 * max(abs(maa$coef))), 
  c("blue", "white", "red")
)
hm1 <- Heatmap((as.matrix(path[,"coef",drop=F])), 
               name = "*FDR<0.05\n**FDR<0.01\nMAA coef",
               cluster_rows = FALSE,
               cluster_columns = T,
               row_names_side = "left",
               show_column_dend = F,
               row_gap = unit(3, "mm"), 
               height = nrow((path[,"coef",drop=F])) * unit(2, "mm"),
               width = ncol((path[,"coef",drop=F])) * unit(2, "mm"),
               rect_gp = gpar(col = "gray", lwd = 0.5),
               row_split = factor(path$Classification,levels=unique(path$Classification)),
               column_names_gp = gpar(fontsize = 5, fontface = "bold"),
               row_names_gp = gpar(fontsize = 5, fontface = "bold"),
               col = col_rho_2,
               column_title_gp = gpar(fontsize = 15),
               row_title_gp = gpar(fontsize = 5, fontface = "bold"),
               row_title_rot = 0,  
               heatmap_legend_param = list(
                 title_gp = gpar(fontsize = 10, fontface = "bold"), 
                 labels_gp = gpar(fontsize = 10),
                 legend_direction = "horizontal", 
                 at = c(-0.1, -0.05, 0, 0.05, 0.1)
               ),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if((path[,"qval",drop=F])[i, j] < 0.05) {
                   if((path[,"qval",drop=F])[i, j] < 0.01) {
                     grid.text("**", x, y, vjust = 0.8, gp = gpar(fontsize = 5))
                   } else {
                     grid.text("*", x, y, vjust = 0.8, gp = gpar(fontsize = 5))
                   }
                 }
               })
set.seed(123)
cairo_pdf(file = paste0( "Fig4C.pdf"),  family="Helvetica",width = 10,height = 10)
ht5 = draw( hm1+hm2, padding = unit(c(80, 50, 2, 2), "mm"), merge_legend=TRUE,main="*FDR<0.05\n**FDR<0.01\nMAA coef"
)
dev.off()

