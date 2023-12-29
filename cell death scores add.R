# add scores of cells


neg_apoptotic<-read.table("D:\\R workspace\\20210812 MDA5 DM new\\calculate_score\\negative_apoptotic_gobp.TXT",
                header = F,sep = "\t")
neg_autophagy<-read.table("D:\\R workspace\\20210812 MDA5 DM new\\calculate_score\\negative_autophagy_gobp.TXT",
                 header = F,sep = "\t")
pos_apoptotic<-read.table("D:\\R workspace\\20210812 MDA5 DM new\\calculate_score\\positive_apoptotic_gobp.TXT",
                          header = F,sep = "\t")
pos_autophagy<-read.table("D:\\R workspace\\20210812 MDA5 DM new\\calculate_score\\positive_autophagy_gobp.TXT",
                          header = F,sep = "\t")
pyroptosis<-read.table("D:\\R workspace\\20210812 MDA5 DM new\\calculate_score\\pyroptosis-reactome.txt",
                          header = F,sep = "\t")

article_apoptosis<-read.table("D:\\R workspace\\20210812 MDA5 DM new\\calculate_score\\apoptosis_article.txt",
                       header = F,sep = "\t")
article_autophagy<-read.table("D:\\R workspace\\20210812 MDA5 DM new\\calculate_score\\autophagy_article.txt",
                              header = F,sep = "\t")
article_necrosis<-read.table("D:\\R workspace\\20210812 MDA5 DM new\\calculate_score\\necrosis_article.txt",
                              header = F,sep = "\t")

neg_apoptotic<-list(intersect(neg_apoptotic[,1],rownames(subT)))
neg_autophagy<-list(intersect(neg_autophagy[,1],rownames(subT)))
pos_apoptotic<-list(intersect(pos_apoptotic[,1],rownames(subT)))
pos_autophagy<-list(intersect(pos_autophagy[,1],rownames(subT)))
pyroptosis<-list(intersect(pyroptosis[,1],rownames(subT)))

article_apoptosis<-list(intersect(article_apoptosis[,1],rownames(subT)))
article_autophagy<-list(intersect(article_autophagy[,1],rownames(subT)))
article_necrosis<-list(intersect(article_necrosis[,1],rownames(subT)))
library(Seurat)
subT<-AddModuleScore(
  object = subT,
  features = pos_apoptotic,
  ctrl = 100,seed="66",
  name = 'pos_apoptotic'
)


subT<-AddModuleScore(
  object = subT,
  features = pos_autophagy,
  name = 'pos_autophagy',
  ctrl = 100,
  seed="66"
)

subT<-AddModuleScore(
  object = subT,
  features = pyroptosis,
  name = 'pyroptosis',
  ctrl = 100,
  seed="66"
)

subT<-AddModuleScore(
  object = subT,
  features = article_necrosis,
  name = 'article_necrosis',
  ctrl = 100,
  seed="66"
)

subT<-AddModuleScore(
  object = subT,
  features = article_autophagy,
  name = 'article_autophagy',
  ctrl = 100,
  seed="66"
)

subT<-AddModuleScore(
  object = subT,
  features = article_apoptosis,
  name = 'article_apoptosis',
  ctrl = 100,
  seed="66"
)

ggplot(data=subT@meta.data,aes(y= pos_autophagy1,x=anno_small,fill=orig.ident))+geom_boxplot()+
  scale_fill_manual(values=groupcolor)+
  stat_compare_means(method = "t.test",hide.ns=T,label="p.signif")+theme_pubr()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggplot(data=subT@meta.data,aes(y= pos_apoptotic1,x=anno_small,fill=orig.ident))+geom_boxplot()+
  scale_fill_manual(values=groupcolor)+
  stat_compare_means(method = "t.test",hide.ns=T,label="p.signif")+theme_pubr()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggplot(data=subT@meta.data,aes(y= neg_autophagy1,x=anno_small,fill=orig.ident))+geom_boxplot()+
  scale_fill_manual(values=groupcolor)+
  stat_compare_means(method = "t.test",hide.ns=T,label="p.signif")+theme_pubr()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggplot(data=subT@meta.data,aes(y= neg_apoptotic1,x=anno_small,fill=orig.ident))+geom_boxplot()+
  scale_fill_manual(values=groupcolor)+
  stat_compare_means(method = "t.test",hide.ns=T,label="p.signif")+theme_pubr()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggplot(data=subT@meta.data,aes(y= pyroptosis1,x=anno_small,fill=orig.ident))+geom_boxplot()+
  scale_fill_manual(values=groupcolor)+
  stat_compare_means(hide.ns=T,label="p.signif")+theme_pubr()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+

ggplot(data=subT@meta.data,aes(y= article_apoptosis1,x=anno_small,fill=orig.ident))+geom_boxplot()+
  scale_fill_manual(values=groupcolor)+
  stat_compare_means(hide.ns=T,label="p.signif")+theme_pubr()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
ggplot(data=subT@meta.data,aes(y= article_autophagy1,x=anno_small,fill=orig.ident))+geom_boxplot()+
  scale_fill_manual(values=groupcolor)+
  stat_compare_means(hide.ns=T,label="p.signif")+theme_pubr()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
ggplot(data=subT@meta.data,aes(y= article_necrosis1,x=anno_small,fill=orig.ident))+geom_boxplot()+
  scale_fill_manual(values=groupcolor)+
  stat_compare_means(hide.ns=T,label="p.signif")+theme_pubr()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


pvalue<-read.table("D:/R workspace/20210812 MDA5 DM new/subT data/4 celldeath score pvalue in subT.TXT",
                   header = T,sep = "\t")
rownames(pvalue)<-levels(subcell$anno_small)
library(pheatmap)
pheatmap(pvalue)

bk <- c(seq(0,0.00009,by=0.00001),seq(0.0001,0.0499,by=0.001),seq(0.05,1,by=0.1))
pheatmap(pvalue,cluster_cols = T,
         cluster_rows = F,border_color = "white",display_numbers=T,cellwidth = 30,cellheight = 30,
         number_format="%.1e",number_color = "white",
         color = c(colorRampPalette(colors = c("#C00000","#fd5f00"))(10),
                   colorRampPalette(colors = c("#fd5f00","#FFFBCB"))(50),
                   rep("#DCDCDC",10)),
         breaks = bk)



g<-intersect(neg_apoptotic$V1,pos_apoptotic$V1)

data<-subT@assays$RNA@data[which(rownames(subT@assays$RNA@data) %in% g),]

data.plot<-data.frame(t=rep(1,dim(data)[1]))
for (celltype in levels(subT$anno_small)){
  HC_POS=which(subT$anno_small==celltype & subT$orig.ident=="HC")
  MDA5_POS=which(subT$anno_small==celltype & subT$orig.ident=="MDA5_patient")
  MEAN=apply(data, 1, function(x){
    HC_mean=mean(x[HC_POS])
    MDA5_mean=mean(x[MDA5_POS])
    return(c(HC_mean,MDA5_mean))})
  data.plot<-cbind(data.plot,t(MEAN))
}
data.plot<-data.plot[,-1]
colnames(data.plot)[seq(1, 25, by = 2)]<-paste("HC",levels(subT$anno_small),sep = "_")
colnames(data.plot)[seq(2, 26, by = 2)]<-paste("MDA5",levels(subT$anno_small),sep = "_")
pheatmap(data.plot)
