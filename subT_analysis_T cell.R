subT <- readRDS("D:/R workspace/20210812 MDA5 DM new/subT data/subT.RDS")

MDA5 <- readRDS("D:/R workspace/20210812 MDA5 DM new/MDA5_final_anno.RDS")
HC <- readRDS("D:/R workspace/20210812 MDA5 DM new/PBMC Healthy Control Data/PBMC_HC_GSE158055.rds")
length(intersect(rownames(MDA5),rownames(HC)))

PBMC <- readRDS("D:/R workspace/20210812 MDA5 DM new/PBMC_merge.RDS")
PBMC@assays$RNA@counts<-PBMC@assays$RNA@counts[intersect(rownames(MDA5),rownames(HC)),]
PBMC@assays$RNA@data<-PBMC@assays$RNA@data[intersect(rownames(MDA5),rownames(HC)),]

Tcolor<-c("#B9EDF8","#39BAE8","#1F6ED4","#20457c",
  "#cff09e","#a8dba8","#79bd9a","#20938b",
  "#f17d80","#dedcee","#ffe600","#F5D769","#ED8A3F")

# only select share genes
# or influence subcluster

subT@assays$RNA@counts<-subT@assays$RNA@counts[intersect(rownames(MDA5),rownames(HC)),]
subT@assays$RNA@data<-subT@assays$RNA@data[intersect(rownames(MDA5),rownames(HC)),]

library(Seurat)
subT <- NormalizeData(subT, normalization.method = "LogNormalize", scale.factor = 10000)
subT <- FindVariableFeatures(subT, selection.method = "vst", nfeatures = 2000)
subT <- ScaleData(subT)
subT <- RunPCA(subT)
# run harmony merge data
library(harmony)
subT<-RunHarmony(subT,"sample_merge")
gc()
subT<- RunUMAP(subT,reduction = "harmony", dims=1:20)
subT<- RunTSNE(subT,reduction = "harmony",dims=1:20)
subT <- FindNeighbors(subT, dims = 1:50,reduction = "harmony")
subT <- FindClusters(subT, resolution = 1)
DimPlot(subT,group.by = "orig.ident",reduction = "tsne",split.by = "orig.ident")
DimPlot(subT,group.by = "sample_merge",reduction = "tsne",split.by = "orig.ident")

DimPlot(subT,group.by = "majortype",reduction = "tsne",label = T)
#DimPlot(subT,group.by = "celltype",reduction = "tsne",label = T)
DimPlot(subT,group.by = "population",reduction = "tsne",label = T)
DimPlot(subT,group.by = "subpopulation",reduction = "tsne",label = T)
DimPlot(subT,reduction = "tsne",label = T)

#mk<-FindMarkers(subT,group.by = "orig.ident",ident.1 = "HC",ident.2 = "MDA5_patient")
mk_all<-FindMarkers(PBMC,group.by = "orig.ident",ident.1 = "HC",ident.2 = "MDA5_patient")

mk=FindAllMarkers(subT,only.pos = T,min.pct = 0.3)
saveRDS(mk_all,"DEGs_HC_MDA5_allPBMC.RDS")

table(subT$seurat_clusters,subT$subpopulation)
table(subT$seurat_clusters,subT$celltype)
table(subT$seurat_clusters)


DotPlot(subT,features = c("CCL5","GZMK","CCR7","GZMA","GNLY","ITGB1","CD4","CD8A","ITGA4","FAM13A","CYLD","TRAT1","S100A4","AHNAK","ANXA1","KLF6","PRF1"))

DotPlot(subT,features = c("CD4","CD8A","CCR7","AQP3","CD69","FOXP3","CCR6","CXCR6","CCL5","PRDM1","GZMK","NKG7","NCAM1","KLRF1","KLRC1","KLRD1"))

subT<-subset(subT,seurat_clusters!=15)
# cellmarkers ###############################
c_name<-data.frame(cluster=as.character(0:14),
                   name=c('CD4_Tcm',"CD4_Tnaive","CD4_Tnaive","NK_FCER1G",
                          "CD8_Tem","CD8_Tnaive","NK_KLRC2","CD8_T_cytotoxic",
                          "CD4_T_cytotoxic","Treg",'CD4_Tem',"Gamma delta T","CD8_Tcm","NKT",
                          'CD4_Tem'))
rownames(c_name)<-c_name$cluster
subT$anno_small<-c_name[subT$seurat_clusters,2]

subT$anno_small<-factor(subT$anno_small,levels = c("CD4_Tnaive",'CD4_Tcm','CD4_Tem',"CD4_T_cytotoxic",
                                                   "CD8_Tnaive",'CD8_Tcm','CD8_Tem',"CD8_T_cytotoxic",
                                                   "Gamma delta T","Treg","NKT","NK_KLRC2","NK_FCER1G"))
library(ggplot2)
library(ggsci)
library(ggforce)
library(ggpubr)  

DimPlot(subT,group.by = "anno_small",reduction = "tsne",label = T)+
     #scale_color_manual(values = pal_simpsons()(16))
     scale_color_manual()+theme_pubr()
  #add circle on tSNE plot
  
  # stat_ellipse(aes(fill="anno_small"),
               # geom = "polygon",
               # alpha=1/5)

data.plot<-subT@reductions$tsne@cell.embeddings
data.plot<-as.data.frame(data.plot)
data.plot$cluster<-subT$anno_small
ggplot(data.plot,mapping = aes(x=tSNE_1,
                         y=tSNE_2,
                         fill=cluster,color=cluster))+theme_pubr()+
  #geom_point(size=1)+scale_color_manual(values = Tcolor)
  stat_ellipse(geom = "polygon",alpha=0.3,type = "t",colour=NA)+ # type: t, norm, euclid
  scale_color_manual(values = Tcolor)+scale_fill_manual(values = Tcolor)


as.data.frame(subT@meta.data) %>%
  group_by(sample_merge) %>%
  summarise(n=n()) ->week.cell

as.data.frame(subT@meta.data) %>%
  group_by(sample_merge,anno_small) %>%
  summarise(cls_n=n()) %>%
  inner_join(week.cell,by=c("sample_merge"="sample_merge")) %>%
  dplyr::mutate(ratio = cls_n/n) ->data.plot

ggplot(data.plot,aes(sample_merge,cls_n,fill=anno_small))+
  geom_bar(stat="identity",position="fill")+theme_pubr()+
  theme(legend.position = "right",legend.title = element_blank())+
  scale_fill_manual(values=c(Tcolor))+
  xlab("")+
  ylab("Cell Compartment")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+


ggplot(data.plot,aes(sample_merge,cls_n,fill=anno_small))+
  geom_bar(stat="identity",position="stack")+theme_pubr()+
  theme(legend.position = "right",legend.title = element_blank())+
  scale_fill_manual(values=Tcolor)+
  xlab("")+
  ylab("Cell Compartment")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))

# stat of diffrent samples
group<-sapply(data.plot$sample_merge,FUN = function(x){
  if(x %in% c("S-HC003","S-HC008","S-HC010")){a="HC"}
  else {a="MDA5_patient"}
  a
})
data.plot$group<-group

p <- ggboxplot(data.plot, x = 'group', y = 'ratio',
               color = 'group', add = "jitter")+
  facet_wrap(~anno_small,nrow = 4)+scale_color_manual(values =groupcolor)
#  Add p-value
p + stat_compare_means(label.y=0.2,label.x = 1)



# density plot####################
library(ggplot2)
library(ggpointdensity)
library(viridis)
data = as.data.frame(subT@reductions$tsne@cell.embeddings)
data$orig.ident<-subT$orig.ident
ggplot(data, mapping = aes(x = tSNE_1, y = tSNE_2)) + facet_wrap(~orig.ident,nrow = 1)+
  geom_pointdensity(size=0.5) + theme_pubr()+scale_color_viridis()

# density plot####################

#  cytokine score and inflammatory score
ifm<-read.table("D:\\R workspace\\20210812 MDA5 DM new\\calculate_score\\Inflammatory score.txt",
                header = F,sep = "\t")
cytk<-read.table("D:\\R workspace\\20210812 MDA5 DM new\\calculate_score\\cytokine score.txt",
                 header = F,sep = "\t")
ifm<-list(intersect(ifm[,1],rownames(subT)))
cytk<-list(intersect(cytk[,1],rownames(subT)))
# 200 inflammatory genes

subT<-AddModuleScore(
  object = subT,
  features = ifm,
  ctrl = 100,
  name = 'inflammatory_score'
)


subT<-AddModuleScore(
  object = subT,
  features = cytk,
  name = 'cytokine_score',
  ctrl = 100,
  seed="66"
)

ggplot(data=subT@meta.data,aes(y= cytokine_score1,x=anno_small,fill=orig.ident))+geom_boxplot()+
  scale_fill_manual(values=groupcolor)+
  stat_compare_means(method = "t.test",hide.ns=T,label="p.signif")+theme_pubr()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggplot(data=subT@meta.data,aes(y= inflammatory_score1,x=anno_small,fill=orig.ident))+geom_boxplot()+
  scale_fill_manual(values=groupcolor)+
  stat_compare_means(method = "t.test",hide.ns=T,label="p.signif")+theme_pubr()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
library(ggplot2)
FeaturePlot(subT,"cytokine_score1",reduction = "tsne")+
  scale_colour_gradient2(low = "#2c7bb6",high = "#ca0020",mid = "#ffffbf",midpoint = 0)

FeaturePlot(subT,"inflammatory_score1",reduction = "tsne")+
  scale_colour_gradient2(low = "#2c7bb6",high = "#ca0020",mid = "#ffffbf",midpoint = 0.05)


library(ggrepel)
subT@meta.data %>% group_by(orig.ident,anno_small) %>% summarise(mean(inflammatory_score1),mean(cytokine_score1))-> score
colnames(score)<-c("condition","population","ifm_score","cytk_score")
ggplot(score,aes(x=ifm_score,y=cytk_score))+geom_point(aes(color=condition))+theme_pubr()+
  geom_label_repel(aes(ifm_score, cytk_score, color=condition,label=population))+
  scale_color_manual(values=groupcolor)

score$ifm_score_norm<-(score$ifm_score-min(score$ifm_score))/(max(score$ifm_score)-min(score$ifm_score))
score$cytk_score_norm<-(score$cytk_score-min(score$cytk_score))/(max(score$cytk_score)-min(score$cytk_score))


score_pair<-data.frame(ifm_HC=(score$ifm_score_norm[1:13])+0.1,ifm_MDA5=(score$ifm_score_norm[14:26])+0.1,
                       cytk_HC=(score$cytk_score_norm[1:13])+0.1,cytk_MDA5=(score$cytk_score_norm[14:26])+0.1)

score_pair %>% mutate(ifm_FC=ifm_MDA5/ifm_HC, cytk_FC=cytk_MDA5/cytk_HC) ->score_pair
                       
rownames(score_pair)<-score$population[1:13]
library(pheatmap)

bk <- c(seq(0.4,0.99,by=0.001),seq(1,6,by=0.01))
pheatmap(score_pair[order(score_pair$cytk_FC,decreasing = T),5:6],cluster_cols = F,
         cluster_rows = F,border_color = "white",display_numbers=T,cellwidth = 30,cellheight = 30,
         color = c(colorRampPalette(colors = c("#2B4162","white"))(600),
                   colorRampPalette(colors = c("white","#C00000"))(501)),
         breaks = bk)


subcell<-subset(subT,anno_small=="CD8_Tnaive")
mk<-FindMarkers(subcell,ident.1 = "HC",ident.2 = "MDA5_patient",group.by = "orig.ident")
mk$gene<-rownames(mk)

mean_gene<-apply(subcell@assays$RNA@data, 1,FUN = function(x) {mean(x)})
#remove mean_gene<0.1 
mk<-mk[mean_gene[rownames(mk)]>0.1,]
mk<-mk[mk$p_val_adj<0.05,]
up<-mk$gene[mk$avg_log2FC<0]
down<-mk$gene[mk$avg_log2FC>0]

# find DEG by myself

HC_exp<-apply(subcell@assays$RNA@data[,colnames(subcell)[which(subcell$orig.ident=="HC")]], 1,
      FUN = function(x) {sum(x>0)})
MDA5_exp<-apply(subcell@assays$RNA@data[,colnames(subcell)[which(subcell$orig.ident=="MDA5_patient")]], 1,
              FUN = function(x) {sum(x>0)})
selecte_gene<-rownames(subcell)[which(HC_exp>(67/2) & MDA5_exp>(28/2))] 

HC_pos<-which(subcell$orig.ident=="HC")
MDA5_pos<-which(subcell$orig.ident=="MDA5_patient")

stat<-apply(subcell@assays$RNA@data[selecte_gene,],1,FUN=function(x){
  res1=wilcox.test(x[HC_pos], x[MDA5_pos], alternative = "two.sided",exact = FALSE)
  res2=mean(x[MDA5_pos])/mean(x[HC_pos])
  return(c(res1$p.value,res2))
})
stat<-as.data.frame(t(stat))
colnames(stat)<-c("p.value","FC")
stat$p.adj<-p.adjust(stat$p.value,"BH")

stat<-stat[order(stat$p.adj),]
stat$logFC<-log2(stat$FC)
library(EnhancedVolcano)
keyvals <- ifelse(
  stat$logFC < -1 & stat$p.adj<0.05, '#84B1ED',
  ifelse(stat$logFC > 1 & stat$p.adj<0.05, '#EC7357',
         '#a3a1a1'))
keyvals[is.na(keyvals)] <- '#a3a1a1'
names(keyvals)[keyvals == '#EC7357'] <- 'MDA5-High'
names(keyvals)[keyvals == '#a3a1a1'] <- 'Not Sig'
names(keyvals)[keyvals == '#84B1ED'] <- 'MDA5-Low'
stat$gene<-rownames(stat)
sel_g<- rownames(stat)[which(abs(stat$logFC) > 1 & stat$p.adj<0.01)]

sel_g<-intersect(intersect(rownames(stat),article_apoptosis[[1]]),
                 rownames(stat)[which(stat$p.adj<0.05 & abs(stat$logFC)>=1)])
EnhancedVolcano(stat,
                lab = rownames(stat),
                x = 'logFC',
                y = 'p.adj',
                selectLab = sel_g,
                xlab = bquote(~Log[2]~ 'fold change'),
                col=c("#a3a1a1","#84B1ED","#60c5ba","#EC7357"),
                colCustom = keyvals,
                axisLabSize = 13,
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 3.0,
                #labvjust = 0.5,
                #labhjust = 1,
                labCol = 'black',
                labFace = 'italic',
                boxedLabels = TRUE,
                colAlpha = 5/5,
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')
##################################################3
test <- bitr((up), 
             toType="ENTREZID", 
             fromType=c("SYMBOL"),  
             OrgDb="org.Hs.eg.db")
test2 <- bitr((down), 
             toType="ENTREZID", 
             fromType=c("SYMBOL"),  
             OrgDb="org.Hs.eg.db")

library(clusterProfiler)
egoBP <- enrichGO(test2$ENTREZID, OrgDb = "org.Hs.eg.db", ont="BP",
                  readable=TRUE,keyType = "ENTREZID")
## use simplify to remove redundant terms
ego2 <- simplify(egoBP, cutoff=0.7, by="p.adjust",
                 select_fun=min)
ego3 <- filter(ego2, p.adjust < 0.01, Count > 10)

ego4 <- mutate(ego3, richFactor = Count / as.numeric
               (sub("/\\d+", "", BgRatio)))
library(ggplot2)
library(forcats)
library(DOSE)
dotplot(ego4,showCategory=10)
ggplot(ego4, showCategory = 10,
       aes(richFactor,fct_reorder(Description, richFactor))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2","#7e62a3"),
                        trans = "log10",
                        guide=guide_colorbar(reverse=TRUE,
                                             order=1)) +
  scale_size_continuous(range=c(1, 8)) +
  theme_dose(12) +
  xlab("Rich Factor") +
  ylab(NULL) +
  ggtitle("GO-BP")+
  theme(axis.text.x= element_text(size=10,angle = 45,hjust = 0.9))+
  theme(legend.key.size = unit(10, "pt"))
# pdf size: 5.4;3.6

library(enrichplot)
heatplot(ego3,showCategory = 4,foldChange = sortFC)+
  scale_fill_gradient(low="#ffda8e",high="#ff5f2e")+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_text(size=4))


cp = list(a.gene=test$ENTREZID, b.gene=test2$ENTREZID)  

xx <- compareCluster(cp, fun="groupGO") 

dotplot(xx,showCategory=10)



##  注释错误
n1=which(subT$anno_small=='CD4_Tem')
n2=which(subT$anno_small=='CD4_T_cytotoxic')

subT$anno_small[n1]='CD4_T_cytotoxic'
subT$anno_small[n2]='CD4_Tem'

VlnPlot(subT,features = c('cytotoxic_score1'),group.by = 'anno_small',pt.size = 0)
