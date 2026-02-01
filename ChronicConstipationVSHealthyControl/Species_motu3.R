library(Maaslin2)

meta5137=as.data.frame(read.csv("metadata.csv"))


fdata=as.data.frame(read.csv("motusp.csv"))
new_old_name=as.data.frame(read.csv("SP_new_old.csv",row.names = 1))
rownames(new_old_name)=new_old_name$new





set.seed(123)
data=fdata[rownames(meta5137),]

identical(rownames(data),rownames(meta5137))
data=data[,(colSums(data>0)/nrow(data))>0.025]
data=data[,-1] #remove unass sp
pr=NULL


diskey=colnames(meta5137)[c(3)]


fixkey=c("Age","Sex","BMI","Vegetable","Fruit","Yogurt", "Bread","Alcohol","Moderate.active")

for (i in diskey){
  set.seed(123)
  
  tempmeta=meta5137
  tempp=NULL
  tempmeta=tempmeta[,c(i,fixkey)]
  tempmeta=na.omit(tempmeta)
  colnames(tempmeta)[1]="test"
  temp=Maaslin2(data[rownames(tempmeta),],tempmeta,
                standardize = T,
                "output",
                fixed_effects = c("test",fixkey),
                min_abundance = 0.0,
                min_prevalence = 0,
                normalization = "NONE",
                plot_heatmap = FALSE,
                plot_scatter = FALSE
  )
  tempp=temp$results
  tempp=tempp[tempp$metadata %in% c("test"),]
  tempp$metadata=gsub("test",i,tempp$metadata)
  tempp$fix=paste(fixkey,collapse = ",")
  colnames(tempp)[1]="feature"
  tempp$old=new_old_name[tempp$feature,]$old
  pr=rbind(pr,tempp)
}

pr=pr[,c(-3,-5,-7)]
write.csv(pr,"./tableS2_motusp.csv",row.names = F)

