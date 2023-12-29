monomac<-subset(MDA5,subset = myanno_new %in% c("mono_Blood","mono_BALF",'macro_m1','macro_m2','macro_m3'))
DimPlot(monomac,group.by = "myanno_new",reduction = "tsne")
monomac$myanno_new<-factor(monomac$myanno_new,levels = c("mono_Blood","mono_BALF",'macro_m1','macro_m2','macro_m3'))
# Destiny diffusionmap ###########################################
library(destiny)
library(Biobase)
library(dplyr)
library(ggsci)
library(ggplot2)
library(gridExtra)
library(ggpubr)

# The input data for "destiny" is ExpressionSet object constructed by Biobase package
ct <-GetAssayData(object = monomac)
monomac <- FindVariableFeatures(monomac, selection.method = "vst", nfeatures = 2000)
ct<-ct[VariableFeatures(monomac),]
ct <- as.ExpressionSet(as.data.frame(t(ct)))


ct$cellType<-factor(monomac$myanno_new,levels = c("mono_Blood","mono_BALF",'macro_m1','macro_m2','macro_m3'))
ct$tissue<-monomac$tissue
ct$orig.ident<-monomac$orig.ident


gc()
dm <- DiffusionMap(ct,n_pcs = 50)

qplot(DC1, DC2, data = dm, colour = cellType) +
  scale_color_manual(values=color_palette)+theme_pubr()+
  theme(legend.position = "right",legend.title =NULL)

qplot(DC1, DC2, data = dm, colour = orig.ident) + scale_color_manual(values=c("#ffb3a7","#f47983","#ff4777","#00bc12","#a1afc9","#44cef6","#4b5cc4"))+
  theme_pubr()+
  theme(legend.position = "right",legend.title =NULL)

qplot(DC1, DC2, data = dm, colour = tissue) + #scale_color_manual(values=color_palette)+
  theme_pubr()+
  theme(legend.position = "right",legend.title =NULL)

dpt1 <- DPT(dm,tips=01)
plot(dpt1,dcs = c(4,1))

dpt2 <- DPT(dm,tips=02)
plot(dpt2,dcs = c(1,4))


dpt3 <- DPT(dm,tips=03)
plot(dpt3,dcs = c(1,4))

plot(dpt, root = 2, paths_to = c(1,3), col_by = 'branch',dcs = c(1,4))





# Slingshot ###########################################

library(SingleCellExperiment)
sim <- SingleCellExperiment(assays = List(counts = monomac@assays$RNA@scale.data))
sum(colnames(monomac)==colnames(sim))
colData(sim)$cellType<-factor(monomac$myanno_new,levels = c("mono_Blood","mono_BALF",'macro_m1','macro_m2','macro_m3'))

pca <- monomac@reductions$pca@cell.embeddings
plot(pca)
DimPlot(monomac,group.by = "myanno_new",reduction = "pca",dims = c(2,3))



library(destiny, quietly = TRUE)
# diffusion
dm <- DiffusionMap(t(assays(sim)$count),n_pcs=50)
library(RColorBrewer)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
plot(dm, c(1,2), pch=16,
     bg="white")

plot.DiffusionMap(dm, dims=c(1,2))

qplot(DC1, DC2, data = dm, colour = cellType) +
  scale_color_manual(values=color_palette)+theme_pubr()+
  theme(legend.position = "right",legend.title =NULL)

rd2 <- dm@eigenvectors

#plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

reducedDims(sim) <- SimpleList(PCA = pca, DiffMap = rd2)


library(slingshot)
sim <- slingshot(sim, clusterLabels = 'cellType', reducedDim = 'PCA')

pesu<-slingPseudotime(sim)
pesu<-as.data.frame(pesu)
pesu$final_id<-monomac$myanno_new
pesu$final_id2<-monomac$orig.ident
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggbeeswarm)
color_palette = c(pal_npg()(9))[c(1,2,7,5)]

ggplot(pesu,aes(y=curve1,x=final_id,colour=final_id))+
  geom_quasirandom()+theme_pubr()+
  theme(legend.position = "none")+
  coord_flip()+
  xlab("")+ylab("pseudotime")
