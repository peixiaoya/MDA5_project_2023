# merge

setwd("D:\\R workspace\\20201208 ZS lung PBMC\\data")

library(Seurat)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "GYM-PBMC-10X5/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc1 <- CreateSeuratObject(counts = pbmc.data, project = "GYM_PBMC", min.cells = 3, min.features = 200)

pbmc.data <- Read10X(data.dir = "TXX-PBMC-10X5/filtered_feature_bc_matrix")
pbmc2 <- CreateSeuratObject(counts = pbmc.data, project = "TXX_PBMC", min.cells = 3, min.features = 200)

pbmc.data <- Read10X(data.dir = "XDJ-PBMC-10X5/filtered_feature_bc_matrix")
pbmc3 <- CreateSeuratObject(counts = pbmc.data, project = "XDJ_PBMC", min.cells = 3, min.features = 200)


# Load the BALF dataset
balf.data <- Read10X(data.dir = "GYM-LUNG-10X5/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
balf1 <- CreateSeuratObject(counts = balf.data, project = "GYM_BALF", min.cells = 3, min.features = 200)

balf.data <- Read10X(data.dir = "TXX-LUNG-10X5/filtered_feature_bc_matrix")
balf2 <- CreateSeuratObject(counts = balf.data, project = "TXX_BALF", min.cells = 3, min.features = 200)

balf.data <- Read10X(data.dir = "XDJ-LUNG-10X5/filtered_feature_bc_matrix")
balf3 <- CreateSeuratObject(counts = balf.data, project = "XDJ_BALF", min.cells = 3, min.features = 200)


rm(pbmc.data);rm(balf.data);rm(raw.data)

length(intersect(intersect(intersect(rownames(pbmc1),rownames(pbmc2)),rownames(pbmc3))),rownames(pbmc4))

gc()

MDA5<- merge(x = pbmc1, y = c(pbmc2,pbmc3,pbmc4,balf1,balf2,balf3))

MDA5[["percent.mt"]] <- PercentageFeatureSet(MDA5, pattern = "^MT-")

MDA5$orig.ident<-factor(MDA5$orig.ident,levels = c("GYM_PBMC","XDJ_PBMC","TXX_PBMC","TXX_PBMC_inactive",
                                                   "GYM_BALF","XDJ_BALF","TXX_BALF"))
setwd("D:/R workspace/20210425 MDA5 DM new")
VlnPlot(MDA5, "nFeature_RNA",pt.size = F,group.by = "orig.ident")
VlnPlot(MDA5, "nCount_RNA",pt.size = F,group.by = "orig.ident")
VlnPlot(MDA5, "percent.mt",pt.size = F,group.by = "orig.ident")

# ribosome
rb.genes <- rownames(MDA5)[grep("^RP[SL]",rownames(MDA5))]
C<-GetAssayData(object = MDA5, slot = "counts")
percent.ribo <- colSums(as.matrix(C[rb.genes,]))/Matrix::colSums(C)*100
MDA5 <- AddMetaData(MDA5, percent.ribo, col.name = "percent.ribo")
VlnPlot(MDA5, "percent.ribo",pt.size = F,group.by = "orig.ident")
 
# because PBMC and BALF ncount_nfeature have different distribution
# so need to different filter standard

PBMC<-merge(x = pbmc1, y = c(pbmc2,pbmc3))
BALF<-merge(x = balf1, y = c(balf2,balf3))

PBMC <- subset(PBMC, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA<20000 & nCount_RNA>500)
BALF <- subset(BALF, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA<30000 & nCount_RNA>500)

MDA5<- merge(x=PBMC,y=BALF)
MDA5$orig.ident<-factor(MDA5$orig.ident,levels = c("GYM_PBMC","XDJ_PBMC","TXX_PBMC","TXX_PBMC_inactive",
                                                   "GYM_BALF","XDJ_BALF","TXX_BALF"))

# remove doublet
# CD3+ CD19+ remove
# CD68+ CD19+ remove
library(dplyr)
a<-as.data.frame(t(sapply(colnames(MDA5),FUN = function(x){
  CD19<-MDA5@assays$RNA@data["CD19",x]
  CD3D<-MDA5@assays$RNA@data["CD3D",x]
  CD68<-MDA5@assays$RNA@data["CD68",x]
  c(CD19,CD3D,CD68)
})))
d1<-rownames(filter(a,V1>0 & V2>0))
d2<-rownames(filter(a,V1>0 & V3>0))
doublet<-c(d1,d2)
cells<-setdiff(colnames(MDA5),doublet)
MDA5<-subset(MDA5,cells = cells)
MDA5<-subset(MDA5,nCount_RNA>500 & percent.mt<10)

table(MDA5$orig.ident)

# gender
VlnPlot(MDA5, features = c('XIST'), slot = "counts", log = TRUE,
        pt.size = F)

MDA5 <- NormalizeData(MDA5, normalization.method = "LogNormalize", scale.factor = 10000)
MDA5 <- FindVariableFeatures(MDA5, selection.method = "vst", nfeatures = 2000)
top10<-head(VariableFeatures(MDA5), 10)
plot1 <- VariableFeaturePlot(MDA5)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

MDA5 <- ScaleData(MDA5)

# cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
MDA5 <- CellCycleScoring(MDA5, s.features = s.genes, g2m.features = g2m.genes)

VlnPlot(MDA5, features = c("S.Score","G2M.Score"))

MDA5 <- RunPCA(MDA5)

# run harmony merge data
library(harmony)
MDA5<-RunHarmony(MDA5,"orig.ident")

gc()

MDA5<- RunUMAP(MDA5,reduction="harmony",dims=1:20)
MDA5 <- FindNeighbors(MDA5, reduction="harmony",dims = 1:20)
MDA5 <- FindClusters(MDA5, resolution = 0.6)
MDA5 <- FindClusters(MDA5, resolution = 0.7)
MDA5 <- FindClusters(MDA5, resolution = 0.8)
MDA5 <- FindClusters(MDA5, resolution = 0.9)
MDA5 <- FindClusters(MDA5, resolution = 1)
gc();DimPlot(MDA5)
saveRDS(MDA5,"MDA5.rds")

setwd("D:\\R workspace\\20210425 MDA5 DM new")
DimPlot(MDA5)
VlnPlot(MDA5,"HBB",pt.size = F)
VlnPlot(MDA5,"GYPA",pt.size = F)
FeaturePlot(MDA5,"HBB")
FeaturePlot(MDA5,"GYPA")

sum(MDA5@assays$RNA@data["HBB",]>0) # 5955
sum(MDA5@assays$RNA@data["GYPA",]>0) # 708
sum(MDA5@assays$RNA@data["HBA1",]>0) # 2343
sum(MDA5@assays$RNA@data["HBA2",]>0) # 4333
sum(MDA5@assays$RNA@data["HBB",]>5) # 2122

# GYM_PBMC has high HBB expression; red blood cells pollution
sum(MDA5@assays$RNA@data["HBB",]>5 & MDA5@assays$RNA@data["GYPA",]>0)
# 700
# HBB>5 mainly occures in GYM_PBMC; and HBB>5 = GYPA>0
# so use HBB>5 as filter criteria
MDA5$HBB<-MDA5@assays$RNA@data["HBB",]
MDA5<-subset(MDA5,subset= HBB<5)


# find markers and annotation

DimPlot(MDA5,label = T)
DimPlot(MDA5,label = T,split.by = "orig.ident",ncol = 4)

MDA5$tissue<-sapply(as.character(MDA5$orig.ident),FUN=function(x){
  t=c()
  if(grepl("PBMC",x)){t=c(t,"PBMC")}
  if(grepl("BALF",x)){t=c(t,"BALF")}
  return(t)
})
MDA5$tissue<-factor(MDA5$tissue,levels = c("PBMC","BALF"))

DimPlot(MDA5,label = T,split.by = "tissue")

FeaturePlot(MDA5,"nCount_RNA",split.by = "tissue")

# Findmarker

markers<-FindAllMarkers(MDA5,min.pct = 0.5)

# T cell 
FeaturePlot(MDA5,"CD3D")
FeaturePlot(MDA5,"CD3E")
FeaturePlot(MDA5,"CD3G")
FeaturePlot(MDA5,"CD247")

FeaturePlot(MDA5,"CD4")
FeaturePlot(MDA5,"CD8A")
FeaturePlot(MDA5,"CD8B")

FeaturePlot(MDA5,"GZMA")

FeaturePlot(MDA5,"HLA-DRA")

FeaturePlot(MDA5,"SELL")

FeaturePlot(MDA5,"CCR7")

FeaturePlot(MDA5,"KLRD1")
FeaturePlot(MDA5,"KLRF1")
FeaturePlot(MDA5,"CD45RA")

# B / ASC
FeaturePlot(MDA5,"CD19")
FeaturePlot(MDA5,"MS4A1")
FeaturePlot(MDA5,"MZB1")
FeaturePlot(MDA5,"JCHAIN")
FeaturePlot(MDA5,"IRF4")
FeaturePlot(MDA5,"IGHG1")
FeaturePlot(MDA5,"PRDM1")
FeaturePlot(MDA5,"XBP1")


# epithelium
FeaturePlot(MDA5,"KRT18")
FeaturePlot(MDA5,"TPPP3")


# DC
FeaturePlot(MDA5,"LILRA4")
FeaturePlot(MDA5,"CD1C")

# monocyte
FeaturePlot(MDA5,"CD14")
FeaturePlot(MDA5,"FCGR3A")
