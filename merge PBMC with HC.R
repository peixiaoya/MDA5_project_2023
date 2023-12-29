pbmc <- readRDS("D:/R workspace/20210812 MDA5 DM new/MDA5_PBMC.RDS")

hc_pbmc<-readRDS("D:/R workspace/20210812 MDA5 DM new/PBMC Healthy Control Data/PBMC_HC_GSE158055.rds")
library(Seurat)
hc_pbmc[["percent.mt"]] <- PercentageFeatureSet(hc_pbmc, pattern = "^MT-")
hc_pbmc<-subset(hc_pbmc,nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA<20000 & nCount_RNA>1500 & percent.mt<=10)

PBMC<-merge(pbmc,hc_pbmc)

PBMC$sample_merge<-apply(PBMC@meta.data,1,FUN = function(x){
  if(x[1]=="SeuratProject"){y=x[20]}
  if(x[1]!="SeuratProject"){y=x[1]}
  return(y)
})
PBMC$sample_merge<-factor(PBMC$sample_merge,levels = c("S-HC003","S-HC008", "S-HC010","GYM_PBMC","TXX_PBMC","XDJ_PBMC"))
PBMC$orig.ident<-sapply(PBMC$orig.ident,FUN = function(x){
  if(x=="SeuratProject"){y="HC"}
  if(x!="SeuratProject"){y="MDA5_patient"}
  return(y)
})
PBMC$orig.ident<-factor(PBMC$orig.ident,levels = c("HC","MDA5_patient"))

PBMC <- NormalizeData(PBMC, normalization.method = "LogNormalize", scale.factor = 10000)
PBMC <- FindVariableFeatures(PBMC, selection.method = "vst", nfeatures = 2000)
PBMC <- ScaleData(PBMC)
PBMC <- RunPCA(PBMC)
# run harmony merge data
library(harmony)
PBMC<-RunHarmony(PBMC,"sample_merge")
gc()
PBMC<- RunUMAP(PBMC,reduction = "harmony", dims=1:20)
PBMC<- RunTSNE(PBMC,reduction = "harmony",dims=1:20)
PBMC <- FindNeighbors(PBMC, dims = 1:20,reduction = "harmony")
PBMC <- FindClusters(PBMC, resolution = 1.0)

DimPlot(PBMC,group.by = "orig.ident",reduction = "tsne")
DimPlot(PBMC,group.by = "sample_merge",reduction = "tsne")

DimPlot(PBMC,group.by = "majortype",reduction = "tsne",label = T)
#DimPlot(PBMC,group.by = "celltype",reduction = "tsne",label = T)
DimPlot(PBMC,group.by = "population",reduction = "tsne",label = T)
DimPlot(PBMC,group.by = "subpopulation",reduction = "tsne",label = T)


PBMC$seurat_clusters<-as.character(PBMC$seurat_clusters)

anno<-data.frame(num=as.character(c(0:19,21,22)),anno=c(rep("T",3),rep("Monocyte",3),
                                          "B",rep("T",3),"B","T","Monocyte",
                                          rep("T",3),"DC","Monocyte","B",
                                          "T","T","DC"))
cell<-data.frame(num=PBMC$seurat_clusters)
library(dplyr)
anno_big<-left_join(cell,anno)

PBMC$anno_big<-anno_big$anno

DimPlot(PBMC,group.by = "anno_big",reduction = "tsne",label = T)+
  scale_color_manual(values = celltypecolor)


DimPlot(PBMC,group.by = "anno_big",reduction = "tsne",label = T,split.by = "orig.ident")+
  scale_color_manual(values = celltypecolor)

# density plot####################
library(ggplot2)
library(ggpointdensity)
library(viridis)
ggplot(data = as.data.frame(PBMC@reductions$tsne@cell.embeddings[which(PBMC$orig.ident=="MDA5_patient" ),]), 
       mapping = aes(x = tSNE_1, y = tSNE_2)) +
  geom_pointdensity(size=1) +
  scale_color_viridis() +
  theme_classic()
# density plot####################


sub_T<-subset(PBMC,anno_big=="T")
sub_B<-subset(PBMC,anno_big=="B")
sub_DC<-subset(PBMC,anno_big=="DC")
sub_mono<-subset(PBMC,anno_big=="Monocyte")


mk<-FindMarkers(PBMC,group.by = "orig.ident",ident.1 = "HC",ident.2 = "MDA5_patient",min.pct = 0.1)
mk$gene<-rownames(mk)

up_sg<-rownames(mk)[mk$avg_log2FC<0 & mk$p_val_adj<0.001 & mk$pct.1>0 & mk$pct.2>0]
down_sg<-rownames(mk)[mk$avg_log2FC>0 & mk$p_val_adj<0.001 & mk$pct.1>0 & mk$pct.2>0]

mk<-mk[which(mk$pct.1>0),]
mk<-mk[which(mk$pct.2>0),]


#########3 unuse 
bulk_deg<-read.csv("D:/R workspace/20210425 MDA5 DM new/PBMC merge/WHP_bulk_result/active_inactive_DEG.csv",header = T,
                   row.names = 1)
up_bk<-rownames(bulk_deg)[bulk_deg$log2FoldChange>0 & bulk_deg$padj<0.05]
down_bk<-rownames(bulk_deg)[bulk_deg$log2FoldChange<0 & bulk_deg$padj<0.05]

intersect(up_bk,up_sg)

intersect(down_bk,down_sg)

# RESULT is not good
# try to use DESEQ2 to analysis DEG

# or DEsingle
# failed; it is too slow
results <- DEsingle(counts = PBMC@assays$RNA@counts, group = PBMC$orig.ident,parallel = TRUE)
# Dividing the DE genes into 3 categories at threshold of FDR < 0.05
results.classified <- DEtype(results = results, threshold = 0.05)
# DEtype subdivides the DE genes found by DEsingle into 3 types: DEs, DEa and DEg
# concrete content please check the DEsingle help doc

# Extract DE genes at threshold of FDR < 0.05
results.sig <- results.classified[results.classified$pvalue.adj.FDR < 0.05, ]
# Extract three types of DE genes separately
results.DEs <- results.sig[results.sig$Type == "DEs", ]
results.DEa <- results.sig[results.sig$Type == "DEa", ]
results.DEg <- results.sig[results.sig$Type == "DEg", ]

# DEs represents differential expression status; 
# DEa represents differential expression abundance; 
# DEg represents general differential expression.


FeaturePlot(PBMC,features = c("CD3E","CD4","CD8A",
                              'CD19', 'CD79A', 'MS4A1',
                              'CD68',"CD14",'LYZ', 
                              "CD1C","FCER1A",'CD1E'),reduction = "tsne",ncol = 3,
            cols =   c("lightgrey", "#d7191c"))
  
