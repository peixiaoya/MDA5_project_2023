# Ref article: SARS-CoV-2 infection triggers profibrotic macrophage responses and lung fibrosis

# using 3 dataset macrophage markers then AddModuleScore
setwd('D:\\R workspace\\20220722_MDA5_replot\\BALF\\Macrophage\\remove c13 14 15\\re-plot')

data1<-read.table('Ayaub ref markers of 2 cluster.txt',header = T,sep='\t')

data2<-read.table('Reyfman ref markers of 5 clusters.txt',header = T,sep='\t')

data3<-read.table('Morse ref markers of 4 cluster.txt',header = T,sep='\t')


#  cytokine score and inflammatory score######################################


for (cluster in unique(data1$Cluster)) {

  genes=list(intersect(data1[which(data1$Cluster==cluster),4],rownames(macro)))
  
  macro<-AddModuleScore(
    object = macro,
    features = genes,
    ctrl = 100,
    name = paste0('Ref_Ayaub_',cluster)
  )
  
}

VlnPlot(macro,features='TGFB1',pt.size=0,)

for (cluster in unique(data2$Cluster)) {
  
  genes=list(intersect(data2[which(data2$Cluster==cluster),4],rownames(macro)))
  
  macro<-AddModuleScore(
    object = macro,
    features = genes,
    ctrl = 100,
    name = paste0('Ref_Reyfman_',cluster)
  )
  
}


for (cluster in unique(data3$Cluster)) {
  
  genes=list(intersect(data3[which(data3$Cluster==cluster),4],rownames(macro)))
  
  macro<-AddModuleScore(
    object = macro,
    features = genes,
    ctrl = 100,
    name = paste0('Ref_Morse_',cluster)
  )
  
}



# 渐变色

colors<-c('#20366b',
          '#4574b2',
          '#7da8c6',
          '#d6ecf3',
          '#fdfdbb',
          '#fdeca3',
          '#fecd83',
          '#f88856',
          '#df4d31',
          '#d53026','#782620')


library(ggpubr)

#VlnPlot(macro,features= ('Ref_Ayaub_FABP4..Mp1'),pt.size = 0,cols = macro_color)+
#  theme_pubr()+theme(axis.text.x = element_text(angle = 60, hjust = 1))

#VlnPlot(macro,features= ('Ref_Ayaub_IPFeMp1'),pt.size = 0,cols = macro_color)+
#  theme_pubr()+theme(axis.text.x = element_text(angle = 60, hjust = 1))

data<-macro@meta.data

data$pos<-apply(data,1,FUN = function(x){
  if (x[29]==levels(data$anno)[1]) {a=1};if (x[29]==levels(data$anno)[6]) {a=6}
  if (x[29]==levels(data$anno)[2]) {a=2};if (x[29]==levels(data$anno)[7]) {a=7}
  if (x[29]==levels(data$anno)[3]) {a=3};if (x[29]==levels(data$anno)[8]) {a=8}
  if (x[29]==levels(data$anno)[4]) {a=4};if (x[29]==levels(data$anno)[9]) {a=9}
  if (x[29]==levels(data$anno)[5]) {a=5};a
})

ggplot(data)+geom_jitter(aes(x=pos,y=Ref_Morse_1.SPP1..Mφ.IPF1,color=Ref_Morse_1.SPP1..Mφ.IPF1),width = 0.33)+
  scale_x_discrete(limits=c(levels(data$anno)))+
  scale_colour_gradientn(colors = colors)+
  geom_violin(aes(x=pos,y=Ref_Morse_1.SPP1..Mφ.IPF1,fill=anno))+
  scale_fill_manual(values = macro_color)+
  theme_classic()+theme(axis.text.x = element_text(angle = 60, hjust = 1))

FeaturePlot(macro,features= ('Ref_Morse_1.SPP1..Mφ.IPF1'))+
  scale_colour_gradientn(colors = colors)


ggplot(data)+geom_jitter(aes(x=pos,y=Ref_Morse_0.FABP4..Mφ1,color=Ref_Morse_0.FABP4..Mφ1),width = 0.33)+
  scale_x_discrete(limits=c(levels(data$anno)))+
  scale_colour_gradientn(colors = colors)+
  geom_violin(aes(x=pos,y=Ref_Morse_0.FABP4..Mφ1,fill=anno))+
  scale_fill_manual(values = macro_color)+
  theme_classic()+theme(axis.text.x = element_text(angle = 60, hjust = 1))

FeaturePlot(macro,features= ('Ref_Morse_0.FABP4..Mφ1'))+
  scale_colour_gradientn(colors = colors)

########################################################################

colnames(macro@meta.data)

# "Ref_Ayaub_IPFeMp1"        "Ref_Ayaub_FABP4..Mp1"

# "Ref_Reyfman_Mφ.Cluster.01" "Ref_Reyfman_Mφ.Cluster.11"
# "Ref_Reyfman_Mφ.Cluster.21" "Ref_Reyfman_Mφ.Cluster.31" "Ref_Reyfman_Mφ.fibrosis1" 

# "Ref_Morse_6.FCN1..Mφ1"     "Ref_Morse_1.SPP1..Mφ1"     "Ref_Morse_1.SPP1..Mφ.IPF1"
# "Ref_Morse_0.FABP4..Mφ1"   


library(MySeuratWrappers) # 可在VlnPlot中使用stacked=T参数
VlnPlot(macro,features = c('TGFB1', 'TGFBI', 'LGMN','CCL18'),pt.size=0,stacked=T,
        cols = macro_color)#+geom_boxplot(width=0.1,fill='white',color='black')

