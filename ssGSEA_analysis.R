# GSEA analysis
library(escape)
library(dittoSeq)

subcell<-subT
# Getting Gene Sets
dataset<-msigdbr::msigdbr_collections()

GS.hallmark <- getGeneSets(species = "Homo sapiens",library = "C5",
                           gene.sets = c("GOBP_NEGATIVE_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY",
                                         "GOBP_POSITIVE_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY",
                                         "GOBP_NEGATIVE_REGULATION_OF_AUTOPHAGY",
                                         "GOBP_POSITIVE_REGULATION_OF_AUTOPHAGY"))

GS.hallmark <- getGeneSets(species = "Homo sapiens",library = "C2")

# Enrichment
res <- enrichIt(obj = subcell, gene.sets = GS.hallmark, groups = 1000, cores = 2)


# add these results back to our Seurat object
subcell <- Seurat::AddMetaData(subcell,res)

subcell$apoptotic<-subcell$GOBP_POSITIVE_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY-subcell$GOBP_NEGATIVE_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY
subcell$autophagy<-subcell$GOBP_POSITIVE_REGULATION_OF_AUTOPHAGY-subcell$GOBP_NEGATIVE_REGULATION_OF_AUTOPHAGY


# Visualizations 
# 1. heatmap
library(dittoSeq)
colors <- colorRampPalette(c("#0348A6", "#7AC5FF", "#C6FDEC", "#FFB433", "#FF4B20"))
dittoHeatmap(subcell, genes = NULL, metas = names(res), 
             annot.by = c("orig.ident"), 
             fontsize = 7, 
             cluster_cols = F,
             heatmap.colors = colors(50))

# 2. violin plot
multi_dittoPlot(subcell, vars = names(res), #"HALLMARK_APOPTOSIS",
                group.by = "anno_small", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Enrichment Scores", 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

# 3. Split Violin Plots
ES2 <- data.frame(subcell[[]], subcell$anno_small)
splitEnrichment(ES2, split = "orig.ident",x.axis = "anno_small",gene.set = "autophagy")

# significance
output <- getSignificance(ES2[,c(1,37:86)], group = "orig.ident", fit = "linear.model")
