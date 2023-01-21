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

