DEG <- readRDS("D:/R workspace/20210812 MDA5 DM new/subT data/CD4_Tcm_gene_expr_stat.RDS")
test <- bitr(rownames(DEG), 
             toType="ENTREZID", 
             fromType="SYMBOL",  
             OrgDb="org.Hs.eg.db")
DEG<-DEG[test$SYMBOL,]
rownames(DEG)<-test$ENTREZID
FC<-DEG$logFC
names(FC)<-rownames(DEG)
sortFC<-sort(FC,decreasing = T)

#######################  GSEA GO ############################
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(forcats)
library(enrichplot)
library(dplyr)
library(ggplot2)
gseBP<-gseGO(sortFC,ont = "BP",keyType = "ENTREZID",pvalueCutoff = 0.05,OrgDb = org.Hs.eg.db)
#gseBP<-na.omit(gseBP)
gseBP2 <- arrange(gseBP, desc(abs(NES))) %>%
  group_by(sign(NES)) %>%
  slice(1:20)

ggplot(gseBP2, showCategory=40,
       aes(NES, fct_reorder(Description, NES),
           fill=qvalues)) +
  geom_col() +
  scale_fill_gradientn(colours=c("#b3eebe","#46bac2", "#371ea3"),
                       guide=guide_colorbar(reverse=TRUE)) +
  theme_dose(12) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("GO-BP")

#######################  GSEA KEGG ############################

kk <- gseKEGG(sortFC, organism = "hsa")
kk2 <- arrange(kk, desc(abs(NES))) %>%
  group_by(sign(NES)) %>%
  slice(1:10)
ggplot(kk2, showCategory=20,
       aes(NES, fct_reorder(Description, NES),
           fill=qvalues)) +
  geom_col() +
  scale_fill_gradientn(colours=c("#b3eebe","#46bac2", "#371ea3"),
                       guide=guide_colorbar(reverse=TRUE)) +
  theme_dose(12) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("KEGG")

gseaplot2(kk,c(9,11,13))

#######################  GSEA WIKIpathway ############################
## downloaded from https://wikipathways-data.wmcloud.org/current/gmt/
gmt <- 'D:\\R workspace\\wipipathway_files\\wikipathways-20210810-gmt-Homo_sapiens.gmt'
wp <- read.gmt.wp(gmt)
ewp <- GSEA(sortFC, TERM2GENE=wp[,c("wpid","gene")], TERM2NAME=wp[,c("wpid", "name")])
ewp2 <- arrange(ewp, desc(abs(NES))) %>%
  group_by(sign(NES)) %>%
  slice(1:20)

# wptest<-gseWP(sort2,organism="mmu")

ggplot(ewp2, showCategory=40,
       aes(NES, fct_reorder(Description, NES),
           fill=qvalues)) +
  geom_col() +
  scale_fill_gradientn(colours=c("#b3eebe","#46bac2", "#371ea3"),
                       guide=guide_colorbar(reverse=TRUE)) +
  theme_dose(12) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("WikiPathways")

gseaplot2(ewp,c(4,7))

ridgeplot(ewp, 10)+scale_fill_gradient(high="#a1bad0",low="#017890")


# #######################  GSEA Reactome Pathway ############################
library(ReactomePA)
Go_Reactomeresult <- gsePathway(sortFC,organism = "human")

gseaplot2(Go_Reactomeresult, 1)
ridgeplot(Go_Reactomeresult, 10)

re2 <- arrange(Go_Reactomeresult, desc(abs(NES))) %>%
  group_by(sign(NES)) %>%
  slice(1:20)

ggplot(re2, showCategory=20,
       aes(NES, fct_reorder(Description, NES),
           fill=qvalues)) +
  geom_col() +
  scale_fill_gradientn(colours=c("#b3eebe","#46bac2", "#371ea3"),
                       guide=guide_colorbar(reverse=TRUE)) +
  theme_dose(12) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("Reactome")

subcell<-subset(subT,anno_small=="CD4_Tnaive")
#gene<-c("CASP1","CASP3","CASP4","CASP5","CASP8","GSDMB","GSMDC","GSDMD","GSDME","GZMA","GZMB")
gene_REACTOME_PYROPTOSIS<-c("BAK1","BAX","CASP1","CASP3","CASP4","CASP5","CHMP2A","CHMP2B","CHMP3","CHMP4A","CHMP4B","CHMP4C","CHMP6","CHMP7","CYCS","ELANE","GSDMD","GSDME","GZMB","HMGB1","IL18","IL1A","IL1B","IRF1","IRF2","TP53","TP63")
DotPlot(subcell,features=gene_REACTOME_PYROPTOSIS,group.by = "orig.ident")
