library(edgeR)
library(limma)
library(BiocManager)
library(devtools)
library(Cairo)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(magick)
#定义颜色
col5 = colorRamp2(c(-2,0,2,6), c("royalblue", "white", "red","black"))
#导入数据
z <- read.csv("D:/Repeats_TE/transcriptid_redoTE/TEdata_output/z.txt",header = T,sep = " ")
SINE <- read.csv("D:/Repeats_TE/transcriptid_redoTE/TEdata_output/SINE.txt",header = T,sep = " ")
LINE <- read.csv("D:/Repeats_TE/transcriptid_redoTE/TEdata_output/LINE.txt",header = T,sep = " ")
LTR <- read.csv("D:/Repeats_TE/transcriptid_redoTE/TEdata_output/LTR.txt",header = T,sep = " ")

NT_KO_vs_WT_up <- read.csv("D:/Repeats_TE/transcriptid_redoTE/TEdata_output/NT_KO_vs_WT_up.txt",header = T,sep = " ")
D5_DAC_KO_vs_WT_up <- read.csv("D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D5_DAC_KO_vs_WT_up.txt",header = T,sep = " ")
D5_BMN_KO_vs_WT_up <- read.csv("D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D5_BMN_KO_vs_WT_up.txt",header = T,sep = " ")
D5_DAC_BMN_KO_vs_WT_up <- read.csv("D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D5_DAC_BMN_KO_vs_WT_up.txt",header = T,sep = " ")
D16_DAC_KO_vs_WT_up <- read.csv("D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D16_DAC_KO_vs_WT_up.txt",header = T,sep = " ")
D16_BMN_KO_vs_WT_up <- read.csv("D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D16_BMN_KO_vs_WT_up.txt",header = T,sep = " ")
D16_DAC_BMN_KO_vs_WT_up <- read.csv("D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D16_DAC_BMN_KO_vs_WT_up.txt",header = T,sep = " ")
D26_DAC_KO_vs_WT_up <- read.csv("D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D26_DAC_KO_vs_WT_up.txt",header = T,sep = " ")
D26_BMN_KO_vs_WT_up <- read.csv("D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D26_BMN_KO_vs_WT_up.txt",header = T,sep = " ")
D26_DAC_BMN_KO_vs_WT_up <- read.csv("D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D26_DAC_BMN_KO_vs_WT_up.txt",header = T,sep = " ")

combined_z <- matrix(nr=551005,nc=20)
combined_z[,1] <- rowMeans(z[,c(1,2)])
combined_z[,2] <- rowMeans(z[,c(3,4)])
combined_z[,3] <- rowMeans(z[,c(5,6)])
combined_z[,4] <- rowMeans(z[,c(7,8)])
combined_z[,5] <- rowMeans(z[,c(9,10)])
combined_z[,6] <- rowMeans(z[,c(11,12)])
combined_z[,7] <- rowMeans(z[,c(13,14)])
combined_z[,8] <- rowMeans(z[,c(15,16)])
combined_z[,9] <- rowMeans(z[,c(17,18)])
combined_z[,10] <- rowMeans(z[,c(19,20)])
combined_z[,11] <- rowMeans(z[,c(21,22)])
combined_z[,12] <- rowMeans(z[,c(23,24)])
combined_z[,13] <- rowMeans(z[,c(25,26)])
combined_z[,14] <- rowMeans(z[,c(27,28)])
combined_z[,15] <- rowMeans(z[,c(29,30)])
combined_z[,16] <- rowMeans(z[,c(31,32)])
combined_z[,17] <- rowMeans(z[,c(33,34)])
combined_z[,18] <- rowMeans(z[,c(35,36)])
combined_z[,19] <- rowMeans(z[,c(37,38)])
combined_z[,20] <- rowMeans(z[,c(39,40)])
rownames(combined_z) <- rownames(z)
colnames(combined_z) <- c("WT_NT","Day5_WT_BMN","Day16_WT_BMN","Day26_WT_BMN","Day5_WT_DAC","Day16_WT_DAC","Day26_WT_DAC","Day5_WT_DAC+BMN","Day16_WT_DAC+BMN","Day26_WT_DAC+BMN","KO_NT","Day5_KO_BMN","Day16_KO_BMN","Day26_KO_BMN","Day5_KO_DAC","Day16_KO_DAC","Day26_KO_DAC","Day5_KO_DAC+BMN","Day16_KO_DAC+BMN","Day26_KO_DAC+BMN")

#Use Z-score to build heatmap out of TE.
#KO vs WT
z_TE_up_KO_vs_WT_DAC <- combined_z[row.names(combined_z)%in%row.names(TE_up_KO_vs_WT_DAC),c(1,5,6,7,11,15,16,17)]
Heatmap(z_TE_up_KO_vs_WT_DAC,
        name="Z-score", 
        col = col5, 
        na_col = "grey", 
        color_space = "LAB",
        column_title = "Activated TEs in 786-O 
(KO vs WT, DAC, 2151)",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 10),
        show_row_name = F,
        show_column_names = T, 
        cluster_columns = F,
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 6),
        use_raster = TRUE,
        raster_quality = 2
)
z_TE_up_KO_vs_WT_DAC_excludeNT <- combined_z[row.names(combined_z)%in%row.names(TE_up_KO_vs_WT_DAC_excludeNT),c(1,5,6,7,11,15,16,17)]
Heatmap(z_TE_up_KO_vs_WT_DAC_excludeNT,
        name="Z-score", 
        col = col5, 
        na_col = "grey", 
        color_space = "LAB",
        column_title = "Activated TEs in 786-O 
(KO vs WT, exclude NT, DAC, 2086)",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 10),
        show_row_name = F,
        show_column_names = T, 
        cluster_columns = F,
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 6),
        use_raster = TRUE,
        raster_quality = 2
)

z_TE_up_KO_vs_WT_BMN <- combined_z[row.names(combined_z)%in%row.names(TE_up_KO_vs_WT_BMN),c(1,2,3,4,11,12,13,14)]
Heatmap(z_TE_up_KO_vs_WT_BMN,
        name="Z-score", 
        col = col5, 
        na_col = "grey", 
        color_space = "LAB",
        column_title = "Activated TEs in 786-O 
(KO vs WT, BMN, 1074)",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 10),
        show_row_name = F,
        show_column_names = T, 
        cluster_columns = F,
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 6),
        use_raster = TRUE,
        raster_quality = 2
)
z_TE_up_KO_vs_WT_BMN_excludeNT <- combined_z[row.names(combined_z)%in%row.names(TE_up_KO_vs_WT_BMN_excludeNT),c(1,2,3,4,11,12,13,14)]
Heatmap(z_TE_up_KO_vs_WT_BMN_excludeNT,
        name="Z-score", 
        col = col5, 
        na_col = "grey", 
        color_space = "LAB",
        column_title = "Activated TEs in 786-O 
(KO vs WT, exclude NT, BMN, 970)",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 10),
        show_row_name = F,
        show_column_names = T, 
        cluster_columns = F,
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 6),
        use_raster = TRUE,
        raster_quality = 2
)


z_TE_up_KO_vs_WT_DAC_BMN <- combined_z[row.names(combined_z)%in%row.names(TE_up_KO_vs_WT_DAC_BMN),c(1,8,9,10,11,18,19,20)]
Heatmap(z_TE_up_KO_vs_WT_DAC_BMN,
        name="Z-score", 
        col = col5, 
        na_col = "grey", 
        color_space = "LAB",
        column_title = "Activated TEs in 786-O 
(KO vs WT, DAC+BMN, 1614)",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 10),
        show_row_name = F,
        show_column_names = T, 
        cluster_columns = F,
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 6),
        use_raster = TRUE,
        raster_quality = 2
)
z_TE_up_KO_vs_WT_DAC_BMN_excludeNT <- combined_z[row.names(combined_z)%in%row.names(TE_up_KO_vs_WT_DAC_BMN_excludeNT),c(1,8,9,10,11,18,19,20)]
Heatmap(z_TE_up_KO_vs_WT_DAC_BMN_excludeNT,
        name="Z-score", 
        col = col5, 
        na_col = "grey", 
        color_space = "LAB",
        column_title = "Activated TEs in 786-O 
(KO vs WT, exclude NT, DAC+BMN, 1566)",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 10),
        show_row_name = F,
        show_column_names = T, 
        cluster_columns = F,
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 6),
        use_raster = TRUE,
        raster_quality = 2
)

z_TE_up_KO_vs_WT_ALL <- combined_z[row.names(combined_z)%in%row.names(TE_up_KO_vs_WT_ALL),]
Heatmap(z_TE_up_KO_vs_WT_ALL,
        name="Z-score", 
        col = col5, 
        na_col = "grey", 
        color_space = "LAB",
        column_title = "Activated TEs in 786-O 
(KO vs WT, ALL, 4094)",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 10),
        show_row_name = F,
        show_column_names = T, 
        cluster_columns = F,
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 6),
        use_raster = TRUE,
        raster_quality = 2
)
z_TE_up_KO_vs_WT_ALL_excludeNT <- combined_z[row.names(combined_z)%in%row.names(TE_up_KO_vs_WT_ALL_excludeNT),]
Heatmap(z_TE_up_KO_vs_WT_ALL_excludeNT,
        name="Z-score", 
        col = col5, 
        na_col = "grey", 
        color_space = "LAB",
        column_title = "Activated TEs in 786-O 
(KO vs WT, exclude NT, ALL, 3985)",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 10),
        show_row_name = F,
        show_column_names = T, 
        cluster_columns = F,
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 6),
        use_raster = TRUE,
        raster_quality = 2
)
