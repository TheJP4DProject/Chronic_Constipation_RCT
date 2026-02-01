
library(vegan)

meta5137=as.data.frame(read.csv("metadata.csv"))
meta5137=meta5137[meta5137$Chronic.contipation==1,]
fdata=as.data.frame(read.csv("MO"))


set.seed(123)
data=fdata[rownames(meta5137),]
mincheck=min(data[data!=0])/2
identical(rownames(data),rownames(meta5137))
data=data[,(colSums(data>0)/nrow(data))>0.025]


prune0w=as.data.frame(read.csv("RCT_0wid.csv"))
prune4w=as.data.frame(read.csv("RCT_4wid.csv"))
prune8w=as.data.frame(read.csv("RCT_8wid.csv"))))
rownames(prune0w)=prune0w$ID
rownames(prune4w)=prune4w$ID
rownames(prune8w)=prune8w$ID
prune0w1key=prune0w[prune0w$Prune==1,]$Metaf
prune0w0key=prune0w[prune0w$Prune==0,]$Metaf

prune4w1key=prune4w[prune4w$Prune==1,]$MetaF
prune4w0key=prune4w[prune4w$Prune==0,]$MetaF

prune8w1key=prune8w[prune8w$Prune==1,]$METAF
prune8w0key=prune8w[prune8w$Prune==0,]$METAF

prune0w1key[! prune0w1key %in% rownames(fdata)]
prune0w0key[! prune0w0key %in% rownames(fdata)]
prune0w1key=prune0w1key[prune0w1key %in% rownames(fdata)]
prune0w0key=prune0w0key[prune0w0key %in% rownames(fdata)]


prune4w1key[! prune4w1key %in% rownames(fdata)]
prune4w0key[! prune4w0key %in% rownames(fdata)]
prune4w1key=prune4w1key[prune4w1key %in% rownames(fdata)]
prune4w0key=prune4w0key[prune4w0key %in% rownames(fdata)]

prune8w1key[! prune8w1key %in% rownames(fdata)]
prune8w0key[! prune8w0key %in% rownames(fdata)]
prune8w1key=prune8w1key[prune8w1key %in% rownames(fdata)]
prune8w0key=prune8w0key[prune8w0key %in% rownames(fdata)]


prune0w1keyid=prune0w[(prune0w$Prune==1)&(prune0w$Metaf %in%prune0w1key),]$ID
prune0w0keyid=prune0w[(prune0w$Prune==0)&(prune0w$Metaf %in%prune0w0key),]$ID

prune4w1keyid=prune4w[(prune4w$Prune==1)&(prune4w$MetaF %in%prune4w1key),]$ID
prune4w0keyid=prune4w[(prune4w$Prune==0)&(prune4w$MetaF %in%prune4w0key),]$ID

prune8w1keyid=prune8w[(prune8w$Prune==1)&(prune8w$METAF %in%prune8w1key),]$ID
prune8w0keyid=prune8w[(prune8w$Prune==0)&(prune8w$METAF %in%prune8w0key),]$ID


check04w1=names(table(c(prune0w1keyid,prune4w1keyid)))[table(c(prune0w1keyid,prune4w1keyid))==2]
check04w0=names(table(c(prune0w0keyid,prune4w0keyid)))[table(c(prune0w0keyid,prune4w0keyid))==2]

check08w1=names(table(c(prune0w1keyid,prune8w1keyid)))[table(c(prune0w1keyid,prune8w1keyid))==2]
check08w0=names(table(c(prune0w0keyid,prune8w0keyid)))[table(c(prune0w0keyid,prune8w0keyid))==2]

check48w1=names(table(c(prune4w1keyid,prune8w1keyid)))[table(c(prune4w1keyid,prune8w1keyid))==2]
check48w0=names(table(c(prune4w0keyid,prune8w0keyid)))[table(c(prune4w0keyid,prune8w0keyid))==2]


fdata=fdata[,colnames(data) ]
pr=NULL
for (i in 1:ncol(data)){
  temp04w1=data.frame("0w"=fdata[prune0w[check04w1,]$Metaf,i],
                      "4w"=fdata[prune4w[check04w1,]$MetaF,i])
  temp04w0=data.frame("0w"=fdata[prune0w[check04w0,]$Metaf,i],
                      "4w"=fdata[prune4w[check04w0,]$MetaF,i])
  
  temp08w1=data.frame("0w"=fdata[prune0w[check08w1,]$Metaf,i],
                      "8w"=fdata[prune8w[check08w1,]$METAF,i])
  temp08w0=data.frame("0w"=fdata[prune0w[check08w0,]$Metaf,i],
                      "8w"=fdata[prune8w[check08w0,]$METAF,i])
  
  temp48w1=data.frame("4w"=fdata[prune4w[check48w1,]$MetaF,i],
                      "8w"=fdata[prune8w[check48w1,]$METAF,i])
  temp48w0=data.frame("4w"=fdata[prune4w[check48w0,]$MetaF,i],
                      "8w"=fdata[prune8w[check48w0,]$METAF,i])
  
  
  result <- data.frame(
    new=colnames(data)[i],
    old=new_old_name[colnames(data)[i],]$old,
    # Group0 0w vs 4w  
    Placebo_04w_log2fc = log2((mean(temp04w0[,2])+mincheck) / (mean(temp04w0[,1])+mincheck)),
    Placebo_04w_pval = wilcox.test(temp04w0[,1], temp04w0[,2], paired = TRUE)$p.value,
    Placebo_04w_fdr=NA,
    
    # Group1 0w vs 4w
    Prune_04w_log2fc = log2((mean(temp04w1[,2])+mincheck) / (mean(temp04w1[,1])+mincheck)),
    Prune_04w_pval = wilcox.test(temp04w1[,1], temp04w1[,2], paired = TRUE)$p.value,
    Prune_04w_fdr=NA,
    
    # Group0 0w vs 8w
    Placebo_08w_log2fc = log2((mean(temp08w0[,2])+mincheck) / (mean(temp08w0[,1])+mincheck)),
    Placebo_08w_pval = wilcox.test(temp08w0[,1], temp08w0[,2], paired = TRUE)$p.value,
    Placebo_08w_fdr=NA,
    # Group1 0w vs 8w
    Prune_08w_log2fc = log2((mean(temp08w1[,2])+mincheck) / (mean(temp08w1[,1])+mincheck)),
    Prune_08w_pval = wilcox.test(temp08w1[,1], temp08w1[,2], paired = TRUE)$p.value,
    Prune_08w_fdr=NA,
    
    # Group0 4w vs 8w
    Placebo_48w_log2fc = log2((mean(temp48w0[,2])+mincheck) / (mean(temp48w0[,1])+mincheck)),
    Placebo_48w_pval = wilcox.test(temp48w0[,1], temp48w0[,2], paired = TRUE)$p.value,
    Placebo_48w_fdr=NA,
    
    # Group1 4w vs 8w
    Prune_48w_log2fc = log2((mean(temp48w1[,2])+mincheck) / (mean(temp48w1[,1])+mincheck)),
    Prune_48w_pval = wilcox.test(temp48w1[,1], temp48w1[,2], paired = TRUE)$p.value,
    Prune_48w_fdr=NA
  )
  
  pr=rbind(pr,result)
  
}
pr$Placebo_04w_fdr = p.adjust(pr$Placebo_04w_pval, "fdr")
pr$Prune_04w_fdr = p.adjust(pr$Prune_04w_pval, "fdr")
pr$Placebo_08w_fdr = p.adjust(pr$Placebo_08w_pval, "fdr")
pr$Prune_08w_fdr = p.adjust(pr$Prune_08w_pval, "fdr")
pr$Placebo_48w_fdr = p.adjust(pr$Placebo_48w_pval, "fdr")
pr$Prune_48w_fdr = p.adjust(pr$Prune_48w_pval, "fdr")


write.csv(pr,"tableS8_RCT_MO.csv")

