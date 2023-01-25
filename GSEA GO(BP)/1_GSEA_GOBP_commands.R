###Top DAC,COMB GSEA GO pathways
library(edgeR)
library(limma)
library(ggplot2)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)#颜色colorRamp2的来源
library(magick)
library(reshape2)
library(data.table)

col1= colorRamp2(c(-2,-1,0,1,2), c("darkgreen","green", "white", "pink","darkred"))

mypathway<- read.csv("C:/GangningLab/RNAseq_EdgeR/GSEA_for_RNAseq_XZ/Top10 DAC,COMB GSEA GO pathways (TvsNT).csv", header = T,sep = ",")
mp_D5_WT_DAC <- mypathway[mypathway$Treatment=="DAC vs. NT Day5" & mypathway$Cell_line=="WT",c(1,2)]
mp_D16_WT_DAC <- mypathway[mypathway$Treatment=="DAC vs. NT Day16" & mypathway$Cell_line=="WT",c(1,2)]
mp_D26_WT_DAC <- mypathway[mypathway$Treatment=="DAC vs. NT Day26" & mypathway$Cell_line=="WT",c(1,2)]
mp_D5_KO_DAC <- mypathway[mypathway$Treatment=="DAC vs. NT Day5" & mypathway$Cell_line=="KO",c(1,2)]
mp_D16_KO_DAC <- mypathway[mypathway$Treatment=="DAC vs. NT Day16" & mypathway$Cell_line=="KO",c(1,2)]
mp_D26_KO_DAC <- mypathway[mypathway$Treatment=="DAC vs. NT Day26" & mypathway$Cell_line=="KO",c(1,2)]
mp_D5_WT_DAC_BMN <- mypathway[mypathway$Treatment=="DAC+BMN vs. NT Day5" & mypathway$Cell_line=="WT",c(1,2)]
mp_D16_WT_DAC_BMN <- mypathway[mypathway$Treatment=="DAC+BMN vs. NT Day16" & mypathway$Cell_line=="WT",c(1,2)]
mp_D26_WT_DAC_BMN <- mypathway[mypathway$Treatment=="DAC+BMN vs. NT Day26" & mypathway$Cell_line=="WT",c(1,2)]
mp_D5_KO_DAC_BMN <- mypathway[mypathway$Treatment=="DAC+BMN vs. NT Day5" & mypathway$Cell_line=="KO",c(1,2)]
mp_D16_KO_DAC_BMN <- mypathway[mypathway$Treatment=="DAC+BMN vs. NT Day16" & mypathway$Cell_line=="KO",c(1,2)]
mp_D26_KO_DAC_BMN <- mypathway[mypathway$Treatment=="DAC+BMN vs. NT Day26" & mypathway$Cell_line=="KO",c(1,2)]

mergemp <- merge(mp_D5_WT_DAC,mp_D16_WT_DAC,by.x = "GS", by.y = "GS",all=TRUE)
mergemp <-merge(mergemp,mp_D26_WT_DAC,by.x = "GS", by.y = "GS",all=TRUE)
mergemp <-merge(mergemp,mp_D5_WT_DAC_BMN,by.x = "GS", by.y = "GS",all=TRUE)
mergemp <-merge(mergemp,mp_D16_WT_DAC_BMN,by.x = "GS", by.y = "GS",all=TRUE)
mergemp <-merge(mergemp,mp_D26_WT_DAC_BMN,by.x = "GS", by.y = "GS",all=TRUE)
mergemp <- merge(mergemp,mp_D5_KO_DAC,by.x = "GS", by.y = "GS",all=TRUE)
mergemp <- merge(mergemp,mp_D16_KO_DAC,by.x = "GS", by.y = "GS",all=TRUE)
mergemp <-merge(mergemp,mp_D26_KO_DAC,by.x = "GS", by.y = "GS",all=TRUE)
mergemp <-merge(mergemp,mp_D5_KO_DAC_BMN,by.x = "GS", by.y = "GS",all=TRUE)
mergemp <-merge(mergemp,mp_D16_KO_DAC_BMN,by.x = "GS", by.y = "GS",all=TRUE)
mergemp <-merge(mergemp,mp_D26_KO_DAC_BMN,by.x = "GS", by.y = "GS",all=TRUE)
rownames(mergemp) <- mergemp$GS;mergemp<-mergemp[,c(2:13)]
colnames(mergemp) <- c("WT DAC Day5","WT DAC Day16","WT DAC Day26","WT DAC+BMN Day5","WT DAC+BMN Day16","WT DAC+BMN Day26","KO DAC Day5","KO DAC Day16","KO DAC Day26","KO DAC+BMN Day5","KO DAC+BMN Day16","KO DAC+BMN Day26")
mergemp[is.na(mergemp)] <- 0

Heatmap(mergemp,
        name="NES", 
        col = col1, 
        na_col = "grey", 
        color_space = "LAB",
        column_title = "GSEA GO(BP)",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 10),
        show_row_name = T,
        show_column_names = T, 
        cluster_columns = F,
        cluster_rows = T,
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 6),
        column_names_rot = 45,
        row_names_gp = gpar(fontsize = 5),
        row_names_max_width = unit(10, "cm"),
        row_names_rot = 0,
        use_raster = TRUE,
        raster_quality = 2
)
mergemp_clear <- mergemp[which(abs(rowSums(mergemp))>6),]
Heatmap(mergemp_clear,
        name="NES", 
        show_heatmap_legend = F,
        col = col1, 
        na_col = "grey", 
        color_space = "LAB",
        column_title = "GSEA GO(BP)",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 10),
        show_row_name = T,
        show_column_names = T, 
        cluster_columns = F,
        cluster_rows = T,
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 6),
        column_names_rot = 45,
        row_names_gp = gpar(fontsize = 5),
        row_names_max_width = unit(10, "cm"),
        row_names_rot = 0,
        use_raster = TRUE,
        raster_quality = 2
)
###screen out top 20 DAC,COMB GSEA GO pathways
write.csv(mergemp,"mergemp.csv")
mergempTOP20up<-read.csv("mergempTOP20up.csv",row.names = 1)
mergempTOP20down<-read.csv("mergempTOP20down.csv",row.names = 1)
#CELL,DRUG,DAY
Heatmap(mergempTOP20up,
        name="NES", 
        show_heatmap_legend = F,
        col = col1, 
        na_col = "grey", 
        color_space = "LAB",
        column_title = "GSEA GO(BP) up-regulated TOP20",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 10),
        show_row_name = T,
        show_column_names = T, 
        cluster_columns = F,
        cluster_rows = T,
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 6),
        column_names_rot = 45,
        row_names_gp = gpar(fontsize = 5),
        row_names_max_width = unit(10, "cm"),
        row_names_rot = 0,
        use_raster = TRUE,
        raster_quality = 2
)

Heatmap(mergempTOP20down,
        name="NES", 
        show_heatmap_legend = F,
        col = col1, 
        na_col = "grey", 
        color_space = "LAB",
        column_title = "GSEA GO(BP) down-regulated TOP20",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 10),
        show_row_name = T,
        show_column_names = T, 
        cluster_columns = F,
        cluster_rows = T,
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 6),
        column_names_rot = 45,
        row_names_gp = gpar(fontsize = 5),
        row_names_max_width = unit(10, "cm"),
        row_names_rot = 0,
        use_raster = TRUE,
        raster_quality = 2
)

###Top BMN GSEA GO pathways
mypathway_BMN<- read.csv("C:/GangningLab/RNAseq_EdgeR/GSEA_for_RNAseq_XZ/Top10 BMN GSEA GO pathways (TvsNT).csv", header = T,sep = ",")
mp_D5_WT_BMN <- mypathway_BMN[mypathway_BMN$Treatment=="BMN vs. NT Day5" & mypathway_BMN$Cell_line=="WT",c(1,2)]
mp_D16_WT_BMN <- mypathway_BMN[mypathway_BMN$Treatment=="BMN vs. NT Day16" & mypathway_BMN$Cell_line=="WT",c(1,2)]
mp_D26_WT_BMN <- mypathway_BMN[mypathway_BMN$Treatment=="BMN vs. NT Day26" & mypathway_BMN$Cell_line=="WT",c(1,2)]
mp_D5_KO_BMN <- mypathway_BMN[mypathway_BMN$Treatment=="BMN vs. NT Day5" & mypathway_BMN$Cell_line=="KO",c(1,2)]
mp_D16_KO_BMN <- mypathway_BMN[mypathway_BMN$Treatment=="BMN vs. NT Day16" & mypathway_BMN$Cell_line=="KO",c(1,2)]
mp_D26_KO_BMN <- mypathway_BMN[mypathway_BMN$Treatment=="BMN vs. NT Day26" & mypathway_BMN$Cell_line=="KO",c(1,2)]
mergemp_BMN <- merge(mp_D5_WT_BMN,mp_D16_WT_BMN,by.x = "GS", by.y = "GS",all=TRUE)
mergemp_BMN <-merge(mergemp_BMN,mp_D26_WT_BMN,by.x = "GS", by.y = "GS",all=TRUE)
mergemp_BMN <- merge(mergemp_BMN,mp_D5_KO_BMN,by.x = "GS", by.y = "GS",all=TRUE)
mergemp_BMN <- merge(mergemp_BMN,mp_D16_KO_BMN,by.x = "GS", by.y = "GS",all=TRUE)
mergemp_BMN <-merge(mergemp_BMN,mp_D26_KO_BMN,by.x = "GS", by.y = "GS",all=TRUE)
rownames(mergemp_BMN) <- mergemp_BMN$GS;mergemp_BMN<-mergemp_BMN[,c(2:7)]
colnames(mergemp_BMN) <- c("WT BMN Day5","WT BMN Day16","WT BMN Day26","KO BMN Day5","KO BMN Day16","KO BMN Day26")
mergemp_BMN[is.na(mergemp_BMN)] <- 0
Heatmap(mergemp_BMN,
        name="NES", 
        col = col1, 
        color_space = "LAB",
        column_title = "GSEA GO(BP)",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 10),
        show_row_name = T,
        show_column_names = T, 
        cluster_columns = F,
        cluster_rows = T,
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 6),
        column_names_rot = 45,
        row_names_gp = gpar(fontsize = 5),
        row_names_max_width = unit(10, "cm"),
        row_names_rot = 0,
        use_raster = TRUE,
        raster_quality = 2
)
###screen out top 20 BMN GSEA GO pathways
write.csv(mergemp_BMN,"mergemp_BMN.csv")
mergempTOP20up_BMN<-read.csv("mergempTOP20up_BMN.csv",row.names = 1)
mergempTOP20down_BMN<-read.csv("mergempTOP20down_BMN.csv",row.names = 1)
#CELL,DRUG,DAY
Heatmap(mergempTOP20up_BMN,
        name="NES", 
        show_heatmap_legend = F,
        col = col1, 
        na_col = "grey", 
        color_space = "LAB",
        column_title = "GSEA GO(BP) up-regulated TOP20 (BMN)",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 10),
        show_row_name = T,
        show_column_names = T, 
        cluster_columns = F,
        cluster_rows = T,
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 6),
        column_names_rot = 45,
        row_names_gp = gpar(fontsize = 5),
        row_names_max_width = unit(10, "cm"),
        row_names_rot = 0,
        use_raster = TRUE,
        raster_quality = 2
)

Heatmap(mergempTOP20down_BMN,
        name="NES", 
        show_heatmap_legend = F,
        col = col1, 
        na_col = "grey", 
        color_space = "LAB",
        column_title = "GSEA GO(BP) down-regulated TOP20 (BMN)",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 10),
        show_row_name = T,
        show_column_names = T, 
        cluster_columns = F,
        cluster_rows = T,
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 6),
        column_names_rot = 45,
        row_names_gp = gpar(fontsize = 5),
        row_names_max_width = unit(10, "cm"),
        row_names_rot = 0,
        use_raster = TRUE,
        raster_quality = 2
)

###gseaplot
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(AnnotationHub)
library(topGO)
library(Rgraphviz)
library(ggplot2)
library(DOSE)
library(stringr)
library(dplyr)

EG2Ensembl=toTable(org.Hs.egENSEMBL)
gene<- bitr(EG2Ensembl$ensembl_id, fromType="ENSEMBL", toType=c("SYMBOL","ENTREZID"), OrgDb="org.Hs.eg.db") 

D5_DAC_BMN_KO_vs_WT_input <- read.csv(file="C:/RNAseq_EdgeR/data_output/edgeRt_D5_DAC_BMN_KO_vs_WT.csv", header = T, sep=',')
D5_DAC_BMN_KO_vs_WT_input<-D5_DAC_BMN_KO_vs_WT_input[,c(1,3)]
D5_DAC_BMN_KO_vs_WT_input<-merge(D5_DAC_BMN_KO_vs_WT_input,gene,by.x="X",by.y="ENSEMBL")
D5_DAC_BMN_KO_vs_WT_GSEA <- D5_DAC_BMN_KO_vs_WT_input$logFC
names(D5_DAC_BMN_KO_vs_WT_GSEA) = D5_DAC_BMN_KO_vs_WT_input$ENTREZID
D5_DAC_BMN_KO_vs_WT_GSEA = sort(D5_DAC_BMN_KO_vs_WT_GSEA, decreasing = TRUE)
D5_DAC_BMN_KO_vs_WT_gseGO <- gseGO(D5_DAC_BMN_KO_vs_WT_GSEA, 'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", nPerm = 1000, minGSSize = 5, maxGSSize = 1000, pvalueCutoff=1)
write.csv(D5_DAC_BMN_KO_vs_WT_gseGO,"D5_DAC_BMN_KO_vs_WT_gseGO.csv")
gseaplot2(D5_DAC_BMN_KO_vs_WT_gseGO,c("GO:0006334","GO:0034728","GO:0034723","GO:0031497"),
          pvalue_table = F,
          title = "GSEA Day5 DAC+BMN KO vs. WT",
          base_size = 4,
          rel_heights = c(1, 0.25, 0.5),
          subplots = 1:3,
          ES_geom = "line")

gseaplot2(D5_DAC_BMN_KO_vs_WT_gseGO,9,pvalue_table = T)

D5_DAC_KO_vs_WT_input <- read.csv(file="C:/RNAseq_EdgeR/data_output/edgeRt_D5_DAC_KO_vs_WT.csv", header = T, sep=',')
D5_DAC_KO_vs_WT_input<-D5_DAC_KO_vs_WT_input[,c(1,3)]
D5_DAC_KO_vs_WT_input<-merge(D5_DAC_KO_vs_WT_input,gene,by.x="X",by.y="ENSEMBL")
D5_DAC_KO_vs_WT_GSEA <- D5_DAC_KO_vs_WT_input$logFC
names(D5_DAC_KO_vs_WT_GSEA) = D5_DAC_KO_vs_WT_input$ENTREZID
D5_DAC_KO_vs_WT_GSEA = sort(D5_DAC_KO_vs_WT_GSEA, decreasing = TRUE)
D5_DAC_KO_vs_WT_gseGO <- gseGO(D5_DAC_KO_vs_WT_GSEA, 'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", nPerm = 1000, minGSSize = 5, maxGSSize = 1000, pvalueCutoff=1)
write.csv(D5_DAC_KO_vs_WT_gseGO,"D5_DAC_KO_vs_WT_gseGO.csv")
gseaplot2(D5_DAC_KO_vs_WT_gseGO,c("GO:0006334","GO:0034728","GO:0034723","GO:0031497"),
          pvalue_table = F,
          title = "GSEA Day5 DAC KO vs. WT",
          base_size = 12,
          rel_heights = c(1, 0.25, 0.5),
          subplots = 1:3,
          ES_geom = "line")

D5_BMN_KO_vs_WT_input <- read.csv(file="C:/RNAseq_EdgeR/data_output/edgeRt_D5_BMN_KO_vs_WT.csv", header = T, sep=',')
D5_BMN_KO_vs_WT_input<-D5_BMN_KO_vs_WT_input[,c(1,3)]
D5_BMN_KO_vs_WT_input<-merge(D5_BMN_KO_vs_WT_input,gene,by.x="X",by.y="ENSEMBL")
D5_BMN_KO_vs_WT_GSEA <- D5_BMN_KO_vs_WT_input$logFC
names(D5_BMN_KO_vs_WT_GSEA) = D5_BMN_KO_vs_WT_input$ENTREZID
D5_BMN_KO_vs_WT_GSEA = sort(D5_BMN_KO_vs_WT_GSEA, decreasing = TRUE)
D5_BMN_KO_vs_WT_gseGO <- gseGO(D5_BMN_KO_vs_WT_GSEA, 'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", nPerm = 1000, minGSSize = 5, maxGSSize = 1000, pvalueCutoff=1)
write.csv(D5_BMN_KO_vs_WT_gseGO,"D5_BMN_KO_vs_WT_gseGO.csv")
gseaplot2(D5_BMN_KO_vs_WT_gseGO,c("GO:0006334","GO:0034728","GO:0034723","GO:0031497"),
          pvalue_table = F,
          title = "GSEA Day5 BMN KO vs. WT",
          base_size = 12,
          rel_heights = c(1, 0.25, 0.5),
          subplots = 1:3,
          ES_geom = "line")

###Volcanoplots
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

GOBP_NUCLEOSOME_ASSEMBLY<-read.csv("C:/RNAseq_EdgeR/GSEA_for_RNAseq_XZ/GOBP_NUCLEOSOME_ASSEMBLY.csv")

edgeRt_WT_D5_DAC_vs_NT<-read.csv("C:/实验学习资料/RNAseq_EdgeR/data_output/edgeRt_WT_D5_DAC_vs_NT.csv")
geneid <- bitr(edgeRt_WT_D5_DAC_vs_NT$X,fromType="ENSEMBL", toType=c("SYMBOL","ENTREZID"), OrgDb="org.Hs.eg.db")
edgeRt_WT_D5_DAC_vs_NT<-merge(edgeRt_WT_D5_DAC_vs_NT,geneid,by.x="X",by.y="ENSEMBL")
edgeRt_KO_D5_DAC_vs_NT<-read.csv("C:/实验学习资料/RNAseq_EdgeR/data_output/edgeRt_KO_D5_DAC_vs_NT.csv")
edgeRt_KO_D5_DAC_vs_NT<-merge(edgeRt_KO_D5_DAC_vs_NT,geneid,by.x="X",by.y="ENSEMBL")
edgeRt_WT_D5_BMN_vs_NT<-read.csv("C:/实验学习资料/RNAseq_EdgeR/data_output/edgeRt_WT_D5_BMN_vs_NT.csv")
edgeRt_WT_D5_BMN_vs_NT<-merge(edgeRt_WT_D5_BMN_vs_NT,geneid,by.x="X",by.y="ENSEMBL")
edgeRt_KO_D5_BMN_vs_NT<-read.csv("C:/实验学习资料/RNAseq_EdgeR/data_output/edgeRt_KO_D5_BMN_vs_NT.csv")
edgeRt_KO_D5_BMN_vs_NT<-merge(edgeRt_KO_D5_BMN_vs_NT,geneid,by.x="X",by.y="ENSEMBL")
edgeRt_WT_D5_DAC_BMN_vs_NT<-read.csv("C:/实验学习资料/RNAseq_EdgeR/data_output/edgeRt_WT_D5_DAC_BMN_vs_NT.csv")
edgeRt_WT_D5_DAC_BMN_vs_NT<-merge(edgeRt_WT_D5_DAC_BMN_vs_NT,geneid,by.x="X",by.y="ENSEMBL")
edgeRt_KO_D5_DAC_BMN_vs_NT<-read.csv("C:/实验学习资料/RNAseq_EdgeR/data_output/edgeRt_KO_D5_DAC_BMN_vs_NT.csv")
edgeRt_KO_D5_DAC_BMN_vs_NT<-merge(edgeRt_KO_D5_DAC_BMN_vs_NT,geneid,by.x="X",by.y="ENSEMBL")
edgeRt_WT_D16_DAC_vs_NT<-read.csv("C:/实验学习资料/RNAseq_EdgeR/data_output/edgeRt_WT_D16_DAC_vs_NT.csv")
edgeRt_WT_D16_DAC_vs_NT<-merge(edgeRt_WT_D16_DAC_vs_NT,geneid,by.x="X",by.y="ENSEMBL")
edgeRt_KO_D16_DAC_vs_NT<-read.csv("C:/实验学习资料/RNAseq_EdgeR/data_output/edgeRt_KO_D16_DAC_vs_NT.csv")
edgeRt_KO_D16_DAC_vs_NT<-merge(edgeRt_KO_D16_DAC_vs_NT,geneid,by.x="X",by.y="ENSEMBL")
edgeRt_WT_D16_BMN_vs_NT<-read.csv("C:/实验学习资料/RNAseq_EdgeR/data_output/edgeRt_WT_D16_BMN_vs_NT.csv")
edgeRt_WT_D16_BMN_vs_NT<-merge(edgeRt_WT_D16_BMN_vs_NT,geneid,by.x="X",by.y="ENSEMBL")
edgeRt_KO_D16_BMN_vs_NT<-read.csv("C:/实验学习资料/RNAseq_EdgeR/data_output/edgeRt_KO_D16_BMN_vs_NT.csv")
edgeRt_KO_D16_BMN_vs_NT<-merge(edgeRt_KO_D16_BMN_vs_NT,geneid,by.x="X",by.y="ENSEMBL")
edgeRt_WT_D16_DAC_BMN_vs_NT<-read.csv("C:/实验学习资料/RNAseq_EdgeR/data_output/edgeRt_WT_D16_DAC_BMN_vs_NT.csv")
edgeRt_WT_D16_DAC_BMN_vs_NT<-merge(edgeRt_WT_D16_DAC_BMN_vs_NT,geneid,by.x="X",by.y="ENSEMBL")
edgeRt_KO_D16_DAC_BMN_vs_NT<-read.csv("C:/实验学习资料/RNAseq_EdgeR/data_output/edgeRt_KO_D16_DAC_BMN_vs_NT.csv")
edgeRt_KO_D16_DAC_BMN_vs_NT<-merge(edgeRt_KO_D16_DAC_BMN_vs_NT,geneid,by.x="X",by.y="ENSEMBL")
edgeRt_WT_D26_DAC_vs_NT<-read.csv("C:/实验学习资料/RNAseq_EdgeR/data_output/edgeRt_WT_D26_DAC_vs_NT.csv")
edgeRt_WT_D26_DAC_vs_NT<-merge(edgeRt_WT_D26_DAC_vs_NT,geneid,by.x="X",by.y="ENSEMBL")
edgeRt_KO_D26_DAC_vs_NT<-read.csv("C:/实验学习资料/RNAseq_EdgeR/data_output/edgeRt_KO_D26_DAC_vs_NT.csv")
edgeRt_KO_D26_DAC_vs_NT<-merge(edgeRt_KO_D26_DAC_vs_NT,geneid,by.x="X",by.y="ENSEMBL")
edgeRt_WT_D26_BMN_vs_NT<-read.csv("C:/实验学习资料/RNAseq_EdgeR/data_output/edgeRt_WT_D26_BMN_vs_NT.csv")
edgeRt_WT_D26_BMN_vs_NT<-merge(edgeRt_WT_D26_BMN_vs_NT,geneid,by.x="X",by.y="ENSEMBL")
edgeRt_KO_D26_BMN_vs_NT<-read.csv("C:/实验学习资料/RNAseq_EdgeR/data_output/edgeRt_KO_D26_BMN_vs_NT.csv")
edgeRt_KO_D26_BMN_vs_NT<-merge(edgeRt_KO_D26_BMN_vs_NT,geneid,by.x="X",by.y="ENSEMBL")
edgeRt_WT_D26_DAC_BMN_vs_NT<-read.csv("C:/实验学习资料/RNAseq_EdgeR/data_output/edgeRt_WT_D26_DAC_BMN_vs_NT.csv")
edgeRt_WT_D26_DAC_BMN_vs_NT<-merge(edgeRt_WT_D26_DAC_BMN_vs_NT,geneid,by.x="X",by.y="ENSEMBL")
edgeRt_KO_D26_DAC_BMN_vs_NT<-read.csv("C:/实验学习资料/RNAseq_EdgeR/data_output/edgeRt_KO_D26_DAC_BMN_vs_NT.csv")
edgeRt_KO_D26_DAC_BMN_vs_NT<-merge(edgeRt_KO_D26_DAC_BMN_vs_NT,geneid,by.x="X",by.y="ENSEMBL")

#GOBP_NUCLEOSOME_ASSEMBLY
#D5
N_WT_D5_DAC_vs_NT <- edgeRt_WT_D5_DAC_vs_NT[edgeRt_WT_D5_DAC_vs_NT$ENTREZID%in% GOBP_NUCLEOSOME_ASSEMBLY$GeneId,]
EnhancedVolcano(N_WT_D5_DAC_vs_NT,
                lab = N_WT_D5_DAC_vs_NT$SYMBOL,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-5, 10),
                title = 'WT Day5 DAC vs NT',
                titleLabSize=10,
                axisLabSize=10,
                subtitle=NULL,
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.7,
                col = c('grey','#00DD77','royalblue','#CC0000'),
                border="full",
                borderWidth=1,
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType="longdash",
                cutoffLineCol="#666666"
)+ theme_bw() + theme(panel.grid=element_blank())

N_KO_D5_DAC_vs_NT <- edgeRt_KO_D5_DAC_vs_NT[edgeRt_KO_D5_DAC_vs_NT$ENTREZID%in% GOBP_NUCLEOSOME_ASSEMBLY$GeneId,]
EnhancedVolcano(N_KO_D5_DAC_vs_NT,
                lab = N_KO_D5_DAC_vs_NT$SYMBOL,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-5, 10),
                title = 'KO Day5 DAC vs NT',
                titleLabSize=10,
                axisLabSize=10,
                subtitle=NULL,
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.7,
                col = c('grey','#00DD77','royalblue','#CC0000'),
                border="full",
                borderWidth=1,
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType="longdash",
                cutoffLineCol="#666666"
)+ theme_bw() + theme(panel.grid=element_blank())

N_WT_D5_BMN_vs_NT <- edgeRt_WT_D5_BMN_vs_NT[edgeRt_WT_D5_BMN_vs_NT$ENTREZID%in% GOBP_NUCLEOSOME_ASSEMBLY$GeneId,]
EnhancedVolcano(N_WT_D5_BMN_vs_NT,
                lab = N_WT_D5_BMN_vs_NT$SYMBOL,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-5, 10),
                title = 'WT Day5 BMN vs NT',
                titleLabSize=10,
                axisLabSize=10,
                subtitle=NULL,
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.7,
                col = c('grey','#00DD77','royalblue','#CC0000'),
                border="full",
                borderWidth=1,
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType="longdash",
                cutoffLineCol="#666666"
)+ theme_bw() + theme(panel.grid=element_blank())

N_KO_D5_BMN_vs_NT <- edgeRt_KO_D5_BMN_vs_NT[edgeRt_KO_D5_BMN_vs_NT$ENTREZID%in% GOBP_NUCLEOSOME_ASSEMBLY$GeneId,]
EnhancedVolcano(N_KO_D5_BMN_vs_NT,
                lab = N_KO_D5_BMN_vs_NT$SYMBOL,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-5, 10),
                title = 'KO Day5 BMN vs NT',
                titleLabSize=10,
                axisLabSize=10,
                subtitle=NULL,
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.7,
                col = c('grey','#00DD77','royalblue','#CC0000'),
                border="full",
                borderWidth=1,
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType="longdash",
                cutoffLineCol="#666666"
)+ theme_bw() + theme(panel.grid=element_blank())

N_WT_D5_DAC_BMN_vs_NT <- edgeRt_WT_D5_DAC_BMN_vs_NT[edgeRt_WT_D5_DAC_BMN_vs_NT$ENTREZID%in% GOBP_NUCLEOSOME_ASSEMBLY$GeneId,]
EnhancedVolcano(N_WT_D5_DAC_BMN_vs_NT,
                lab = N_WT_D5_DAC_BMN_vs_NT$SYMBOL,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-5, 10),
                title = 'WT Day5 DAC_BMN vs NT',
                titleLabSize=10,
                axisLabSize=10,
                subtitle=NULL,
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.7,
                col = c('grey','#00DD77','royalblue','#CC0000'),
                border="full",
                borderWidth=1,
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType="longdash",
                cutoffLineCol="#666666"
)+ theme_bw() + theme(panel.grid=element_blank())

N_KO_D5_DAC_BMN_vs_NT <- edgeRt_KO_D5_DAC_BMN_vs_NT[edgeRt_KO_D5_DAC_BMN_vs_NT$ENTREZID%in% GOBP_NUCLEOSOME_ASSEMBLY$GeneId,]
EnhancedVolcano(N_KO_D5_DAC_BMN_vs_NT,
                lab = N_KO_D5_DAC_BMN_vs_NT$SYMBOL,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-5, 10),
                title = 'KO Day5 DAC_BMN vs NT',
                titleLabSize=10,
                axisLabSize=10,
                subtitle=NULL,
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.7,
                col = c('grey','#00DD77','royalblue','#CC0000'),
                border="full",
                borderWidth=1,
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType="longdash",
                cutoffLineCol="#666666"
)+ theme_bw() + theme(panel.grid=element_blank())

#D16
N_WT_D16_DAC_vs_NT <- edgeRt_WT_D16_DAC_vs_NT[edgeRt_WT_D16_DAC_vs_NT$ENTREZID%in% GOBP_NUCLEOSOME_ASSEMBLY$GeneId,]
EnhancedVolcano(N_WT_D16_DAC_vs_NT,
                lab = N_WT_D16_DAC_vs_NT$SYMBOL,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-12.5, 12.5),
                title = 'WT Day16 DAC vs NT',
                titleLabSize=10,
                axisLabSize=10,
                subtitle=NULL,
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.7,
                col = c('grey','#00DD77','royalblue','#CC0000'),
                border="full",
                borderWidth=1,
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType="longdash",
                cutoffLineCol="#666666"
)+ theme_bw() + theme(panel.grid=element_blank())

N_KO_D16_DAC_vs_NT <- edgeRt_KO_D16_DAC_vs_NT[edgeRt_KO_D16_DAC_vs_NT$ENTREZID%in% GOBP_NUCLEOSOME_ASSEMBLY$GeneId,]
EnhancedVolcano(N_KO_D16_DAC_vs_NT,
                lab = N_KO_D16_DAC_vs_NT$SYMBOL,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-12.5, 12.5),
                title = 'KO Day16 DAC vs NT',
                titleLabSize=10,
                axisLabSize=10,
                subtitle=NULL,
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.7,
                col = c('grey','#00DD77','royalblue','#CC0000'),
                border="full",
                borderWidth=1,
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType="longdash",
                cutoffLineCol="#666666"
)+ theme_bw() + theme(panel.grid=element_blank())

N_WT_D16_BMN_vs_NT <- edgeRt_WT_D16_BMN_vs_NT[edgeRt_WT_D16_BMN_vs_NT$ENTREZID%in% GOBP_NUCLEOSOME_ASSEMBLY$GeneId,]
EnhancedVolcano(N_WT_D16_BMN_vs_NT,
                lab = N_WT_D16_BMN_vs_NT$SYMBOL,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-12.5, 12.5),
                title = 'WT Day16 BMN vs NT',
                titleLabSize=10,
                axisLabSize=10,
                subtitle=NULL,
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.7,
                col = c('grey','#00DD77','royalblue','#CC0000'),
                border="full",
                borderWidth=1,
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType="longdash",
                cutoffLineCol="#666666"
)+ theme_bw() + theme(panel.grid=element_blank())

N_KO_D16_BMN_vs_NT <- edgeRt_KO_D16_BMN_vs_NT[edgeRt_KO_D16_BMN_vs_NT$ENTREZID%in% GOBP_NUCLEOSOME_ASSEMBLY$GeneId,]
EnhancedVolcano(N_KO_D16_BMN_vs_NT,
                lab = N_KO_D16_BMN_vs_NT$SYMBOL,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-12.5, 12.5),
                title = 'KO Day16 BMN vs NT',
                titleLabSize=10,
                axisLabSize=10,
                subtitle=NULL,
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.7,
                col = c('grey','#00DD77','royalblue','#CC0000'),
                border="full",
                borderWidth=1,
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType="longdash",
                cutoffLineCol="#666666"
)+ theme_bw() + theme(panel.grid=element_blank())

N_WT_D16_DAC_BMN_vs_NT <- edgeRt_WT_D16_DAC_BMN_vs_NT[edgeRt_WT_D16_DAC_BMN_vs_NT$ENTREZID%in% GOBP_NUCLEOSOME_ASSEMBLY$GeneId,]
EnhancedVolcano(N_WT_D16_DAC_BMN_vs_NT,
                lab = N_WT_D16_DAC_BMN_vs_NT$SYMBOL,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-12.5, 12.5),
                title = 'WT Day16 DAC_BMN vs NT',
                titleLabSize=10,
                axisLabSize=10,
                subtitle=NULL,
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.7,
                col = c('grey','#00DD77','royalblue','#CC0000'),
                border="full",
                borderWidth=1,
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType="longdash",
                cutoffLineCol="#666666"
)+ theme_bw() + theme(panel.grid=element_blank())

N_KO_D16_DAC_BMN_vs_NT <- edgeRt_KO_D16_DAC_BMN_vs_NT[edgeRt_KO_D16_DAC_BMN_vs_NT$ENTREZID%in% GOBP_NUCLEOSOME_ASSEMBLY$GeneId,]
EnhancedVolcano(N_KO_D16_DAC_BMN_vs_NT,
                lab = N_KO_D16_DAC_BMN_vs_NT$SYMBOL,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-12.5, 12.5),
                title = 'KO Day16 DAC_BMN vs NT',
                titleLabSize=10,
                axisLabSize=10,
                subtitle=NULL,
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.7,
                col = c('grey','#00DD77','royalblue','#CC0000'),
                border="full",
                borderWidth=1,
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType="longdash",
                cutoffLineCol="#666666"
)+ theme_bw() + theme(panel.grid=element_blank())

#D26
N_WT_D26_DAC_vs_NT <- edgeRt_WT_D26_DAC_vs_NT[edgeRt_WT_D26_DAC_vs_NT$ENTREZID%in% GOBP_NUCLEOSOME_ASSEMBLY$GeneId,]
EnhancedVolcano(N_WT_D26_DAC_vs_NT,
                lab = N_WT_D26_DAC_vs_NT$SYMBOL,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-12.5, 12.5),
                title = 'WT Day26 DAC vs NT',
                titleLabSize=10,
                axisLabSize=10,
                subtitle=NULL,
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.7,
                col = c('grey','#00DD77','royalblue','#CC0000'),
                border="full",
                borderWidth=1,
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType="longdash",
                cutoffLineCol="#666666"
)+ theme_bw() + theme(panel.grid=element_blank())

N_KO_D26_DAC_vs_NT <- edgeRt_KO_D26_DAC_vs_NT[edgeRt_KO_D26_DAC_vs_NT$ENTREZID%in% GOBP_NUCLEOSOME_ASSEMBLY$GeneId,]
EnhancedVolcano(N_KO_D26_DAC_vs_NT,
                lab = N_KO_D26_DAC_vs_NT$SYMBOL,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-12.5, 12.5),
                title = 'KO Day26 DAC vs NT',
                titleLabSize=10,
                axisLabSize=10,
                subtitle=NULL,
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.7,
                col = c('grey','#00DD77','royalblue','#CC0000'),
                border="full",
                borderWidth=1,
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType="longdash",
                cutoffLineCol="#666666"
)+ theme_bw() + theme(panel.grid=element_blank())

N_WT_D26_BMN_vs_NT <- edgeRt_WT_D26_BMN_vs_NT[edgeRt_WT_D26_BMN_vs_NT$ENTREZID%in% GOBP_NUCLEOSOME_ASSEMBLY$GeneId,]
EnhancedVolcano(N_WT_D26_BMN_vs_NT,
                lab = N_WT_D26_BMN_vs_NT$SYMBOL,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-12.5, 12.5),
                title = 'WT Day26 BMN vs NT',
                titleLabSize=10,
                axisLabSize=10,
                subtitle=NULL,
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.7,
                col = c('grey','#00DD77','royalblue','#CC0000'),
                border="full",
                borderWidth=1,
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType="longdash",
                cutoffLineCol="#666666"
)+ theme_bw() + theme(panel.grid=element_blank())

N_KO_D26_BMN_vs_NT <- edgeRt_KO_D26_BMN_vs_NT[edgeRt_KO_D26_BMN_vs_NT$ENTREZID%in% GOBP_NUCLEOSOME_ASSEMBLY$GeneId,]
EnhancedVolcano(N_KO_D26_BMN_vs_NT,
                lab = N_KO_D26_BMN_vs_NT$SYMBOL,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-12.5, 12.5),
                title = 'KO Day26 BMN vs NT',
                titleLabSize=10,
                axisLabSize=10,
                subtitle=NULL,
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.7,
                col = c('grey','#00DD77','royalblue','#CC0000'),
                border="full",
                borderWidth=1,
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType="longdash",
                cutoffLineCol="#666666"
)+ theme_bw() + theme(panel.grid=element_blank())

N_WT_D26_DAC_BMN_vs_NT <- edgeRt_WT_D26_DAC_BMN_vs_NT[edgeRt_WT_D26_DAC_BMN_vs_NT$ENTREZID%in% GOBP_NUCLEOSOME_ASSEMBLY$GeneId,]
EnhancedVolcano(N_WT_D26_DAC_BMN_vs_NT,
                lab = N_WT_D26_DAC_BMN_vs_NT$SYMBOL,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-12.5, 12.5),
                title = 'WT Day26 DAC_BMN vs NT',
                titleLabSize=10,
                axisLabSize=10,
                subtitle=NULL,
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.7,
                col = c('grey','#00DD77','royalblue','#CC0000'),
                border="full",
                borderWidth=1,
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType="longdash",
                cutoffLineCol="#666666"
)+ theme_bw() + theme(panel.grid=element_blank())

N_KO_D26_DAC_BMN_vs_NT <- edgeRt_KO_D26_DAC_BMN_vs_NT[edgeRt_KO_D26_DAC_BMN_vs_NT$ENTREZID%in% GOBP_NUCLEOSOME_ASSEMBLY$GeneId,]
EnhancedVolcano(N_KO_D26_DAC_BMN_vs_NT,
                lab = N_KO_D26_DAC_BMN_vs_NT$SYMBOL,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-12.5, 12.5),
                title = 'KO Day26 DAC_BMN vs NT',
                titleLabSize=10,
                axisLabSize=10,
                subtitle=NULL,
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 0.7,
                col = c('grey','#00DD77','royalblue','#CC0000'),
                border="full",
                borderWidth=1,
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType="longdash",
                cutoffLineCol="#666666"
)+ theme_bw() + theme(panel.grid=element_blank())
