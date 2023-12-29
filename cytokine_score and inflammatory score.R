#  cytokine score and inflammatory score
ifm<-read.table("D:\\R workspace\\20210812 MDA5 DM new\\calculate_score\\Inflammatory score.txt",
                header = F,sep = "\t")
cytk<-read.table("D:\\R workspace\\20210812 MDA5 DM new\\calculate_score\\cytokine score.txt",
                 header = F,sep = "\t")
ifm<-list(intersect(ifm[,1],rownames(PBMC)))
cytk<-list(intersect(cytk[,1],rownames(PBMC)))
# 200 inflammatory genes

PBMC<-AddModuleScore(
  object = PBMC,
  features = ifm,
  ctrl = 100,
  name = 'inflammatory_score'
)


PBMC<-AddModuleScore(
  object = PBMC,
  features = cytk,
  name = 'cytokine_score',
  ctrl = 100,
  seed="66"
)

ggplot(data=PBMC@meta.data,aes(y= cytokine_score1,x=anno_big,fill=orig.ident))+geom_boxplot()+
         scale_fill_manual(values=samplecolor)+stat_compare_means(method = "t.test")+theme_pubr()
ggplot(data=PBMC@meta.data,aes(y= inflammatory_score1,x=anno_big,fill=orig.ident))+geom_boxplot()+
  scale_fill_manual(values=samplecolor)+stat_compare_means(method = "t.test")+theme_pubr()
library(ggplot2)
FeaturePlot(PBMC,"cytokine_score1",reduction = "tsne")+
  scale_colour_gradient2(low = "#2c7bb6",high = "#ca0020",mid = "#ffffbf",midpoint = 0.7)

FeaturePlot(PBMC,"inflammatory_score1",reduction = "tsne")+
  scale_colour_gradient2(low = "#2c7bb6",high = "#ca0020",mid = "#ffffbf",midpoint = 0.15)


library(ggrepel)
PBMC@meta.data %>% group_by(orig.ident,anno_big) %>% summarise(mean(inflammatory_score1),mean(cytokine_score1))-> score
colnames(score)<-c("condition","population","ifm_score","cytk_score")
ggplot(score,aes(x=ifm_score,y=cytk_score))+geom_point(aes(color=condition))+theme_pubr()+
  geom_label_repel(aes(ifm_score, cytk_score, color=condition,label=population))+
  scale_color_manual(values=samplecolor)
