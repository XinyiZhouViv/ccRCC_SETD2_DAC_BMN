library(edgeR)
library(limma)
library(BiocManager)
library(devtools)
library(Cairo)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(magick)
col5 = colorRamp2(c(-2,0,2,5), c("royalblue", "white", "red","black"))
myTEdata <- read.table("D:/Repeats_TE/transcriptid_redoTE/output_Repeats_transcriptid_classid.txt", header=T)
row.names(myTEdata) <- myTEdata$Geneid
colnames(myTEdata) <- c("Geneid","Classid","Chr","Start","End","Strand","Length","1_WT_NT","1_Day5_KO_BMN","1_Day16_KO_BMN","1_Day26_KO_BMN","1_KO_NT","1_Day5_KO_DAC+BMN","1_Day16_KO_DAC+BMN","1_Day26_KO_DAC+BMN","2_WT_NT","2_Day5_WT_DAC","2_Day16_WT_DAC","2_Day26_WT_DAC","2_KO_NT","2_Day5_WT_BMN","2_Day16_WT_BMN","2_Day26_WT_BMN","1_Day5_WT_DAC","2_Day5_WT_DAC+BMN","1_Day16_WT_DAC","2_Day16_WT_DAC+BMN","1_Day26_WT_DAC","2_Day26_WT_DAC+BMN","1_Day5_WT_BMN","2_Day5_KO_DAC","1_Day16_WT_BMN","2_Day16_KO_DAC","1_Day26_WT_BMN","2_Day26_KO_DAC","1_Day5_WT_DAC+BMN","2_Day5_KO_BMN","1_Day16_WT_DAC+BMN","2_Day16_KO_BMN","1_Day26_WT_DAC+BMN","2_Day26_KO_BMN","1_Day5_KO_DAC","2_Day5_KO_DAC+BMN","1_Day16_KO_DAC","2_Day16_KO_DAC+BMN","1_Day26_KO_DAC","2_Day26_KO_DAC+BMN")
myTEdata <- myTEdata[,c(1:7,8,16,30,21,32,22,34,23,24,17,26,18,28,19,36,25,38,27,40,29,12,20,9,37,10,39,11,41,42,31,44,33,46,35,13,43,14,45,15,47)]
x <- myTEdata[ ,8:47]
group = rep(1:20,each=2)
y <- DGEList(counts=x,group=group, lib.size = colSums(x), genes=data.frame(Length=myTEdata$Length))
keep_cpm <- rowSums(cpm(y)>1) >= 2
y<-y[keep_cpm, ,keep.lib.sizes=FALSE] #gene number from 4693511 to 551005
y <- calcNormFactors(y)
y$samples #check normalization factor
TE_list <- y$genes
#z-score
x <- x[row.names(x)%in%row.names(TE_list),]
x<-transform(x,mean=rowMeans(x))
x <-transform(x,sd=apply(x[,1:40],1,sd))
z<-x
df <- matrix(nrow = 551005, ncol = 40)

for (j in 27:40){
  for (i in 1:551005){
    
    df[i,j] = (z[i,j] - z[i,41] ) / z[i,42]
    
  }
}
z<-df
row.names(z)=row.names(x)
colnames(z)=colnames(x[,1:40])
#SINE, LINE, LTR list
SINE <- myTEdata[myTEdata$Classid=="SINE",c(1,2)]
LINE <- myTEdata[myTEdata$Classid=="LINE",c(1,2)]
LTR <- myTEdata[myTEdata$Classid=="LTR",c(1,2)]

#KO vs WT
NT_KO_vs_WT <- y[,c(1,2,21,22)]
D5_DAC_KO_vs_WT <- y[,c(9,10,29,30)]
D5_BMN_KO_vs_WT <- y[,c(3,4,23,24)]
D5_DAC_BMN_KO_vs_WT <- y[,c(15,16,35,36)]
D16_DAC_KO_vs_WT <- y[,c(11,12,31,32)]
D16_BMN_KO_vs_WT <- y[,c(5,6,25,26)]
D16_DAC_BMN_KO_vs_WT <- y[,c(17,18,37,38)]
D26_DAC_KO_vs_WT <- y[,c(13,14,33,34)]
D26_BMN_KO_vs_WT <- y[,c(7,8,27,28)]
D26_DAC_BMN_KO_vs_WT <- y[,c(19,20,39,40)]

z_NT_KO_vs_WT <- estimateDisp(NT_KO_vs_WT)
print(paste("sqrt of common dispersion (BCV) of","KO_NT", "Vs", "WT_NT", "is", sqrt(z_NT_KO_vs_WT$common.dispersion)))
#[1] "sqrt of common dispersion (BCV) of KO_NT Vs WT_NT is 0.37095353858701"
et_NT_KO_vs_WT <- exactTest(z_NT_KO_vs_WT)
results_edgeR_NT_KO_vs_WT <- topTags(et_NT_KO_vs_WT, n = nrow(x), sort.by = "none")
head(results_edgeR_NT_KO_vs_WT$table)
edgeRt_NT_KO_vs_WT <- results_edgeR_NT_KO_vs_WT$table

z_D5_DAC_KO_vs_WT <- estimateDisp(D5_DAC_KO_vs_WT)
print(paste("sqrt of common dispersion (BCV) of","D5_DAC_KO", "Vs", "D5_DAC_WT", "is", sqrt(z_D5_DAC_KO_vs_WT$common.dispersion)))
#[1] "sqrt of common dispersion (BCV) of D5_DAC_KO Vs D5_DAC_WT is 0.354082181796309"
et_D5_DAC_KO_vs_WT <- exactTest(z_D5_DAC_KO_vs_WT)
results_edgeR_D5_DAC_KO_vs_WT <- topTags(et_D5_DAC_KO_vs_WT, n = nrow(x), sort.by = "none")
head(results_edgeR_D5_DAC_KO_vs_WT$table)
edgeRt_D5_DAC_KO_vs_WT <- results_edgeR_D5_DAC_KO_vs_WT$table

z_D5_BMN_KO_vs_WT <- estimateDisp(D5_BMN_KO_vs_WT)
print(paste("sqrt of common dispersion (BCV) of","D5_BMN_KO", "Vs", "D5_BMN_WT", "is", sqrt(z_D5_BMN_KO_vs_WT$common.dispersion)))
#[1] "sqrt of common dispersion (BCV) of D5_BMN_KO Vs D5_BMN_WT is 0.322993852425662"
et_D5_BMN_KO_vs_WT <- exactTest(z_D5_BMN_KO_vs_WT)
results_edgeR_D5_BMN_KO_vs_WT <- topTags(et_D5_BMN_KO_vs_WT, n = nrow(x), sort.by = "none")
head(results_edgeR_D5_BMN_KO_vs_WT$table)
edgeRt_D5_BMN_KO_vs_WT <- results_edgeR_D5_BMN_KO_vs_WT$table

z_D5_DAC_BMN_KO_vs_WT <- estimateDisp(D5_DAC_BMN_KO_vs_WT)
print(paste("sqrt of common dispersion (BCV) of","D5_DAC_BMN_KO", "Vs", "D5_DAC_BMN_WT", "is", sqrt(z_D5_DAC_BMN_KO_vs_WT$common.dispersion)))
#[1] "sqrt of common dispersion (BCV) of D5_DAC_BMN_KO Vs D5_DAC_BMN_WT is 0.308391727931952"
et_D5_DAC_BMN_KO_vs_WT <- exactTest(z_D5_DAC_BMN_KO_vs_WT)
results_edgeR_D5_DAC_BMN_KO_vs_WT <- topTags(et_D5_DAC_BMN_KO_vs_WT, n = nrow(x), sort.by = "none")
head(results_edgeR_D5_DAC_BMN_KO_vs_WT$table)
edgeRt_D5_DAC_BMN_KO_vs_WT <- results_edgeR_D5_DAC_BMN_KO_vs_WT$table

z_D16_DAC_KO_vs_WT <- estimateDisp(D16_DAC_KO_vs_WT)
print(paste("sqrt of common dispersion (BCV) of","D16_DAC_KO", "Vs", "D16_DAC_WT", "is", sqrt(z_D16_DAC_KO_vs_WT$common.dispersion)))
#[1] "sqrt of common dispersion (BCV) of D16_DAC_KO Vs D16_DAC_WT is 0.307935486921275"
et_D16_DAC_KO_vs_WT <- exactTest(z_D16_DAC_KO_vs_WT)
results_edgeR_D16_DAC_KO_vs_WT <- topTags(et_D16_DAC_KO_vs_WT, n = nrow(x), sort.by = "none")
head(results_edgeR_D16_DAC_KO_vs_WT$table)
edgeRt_D16_DAC_KO_vs_WT <- results_edgeR_D16_DAC_KO_vs_WT$table

z_D16_BMN_KO_vs_WT <- estimateDisp(D16_BMN_KO_vs_WT)
print(paste("sqrt of common dispersion (BCV) of","D16_BMN_KO", "Vs", "D16_BMN_WT", "is", sqrt(z_D16_BMN_KO_vs_WT$common.dispersion)))
#[1] "sqrt of common dispersion (BCV) of D16_BMN_KO Vs D16_BMN_WT is 0.312083739636053"
et_D16_BMN_KO_vs_WT <- exactTest(z_D16_BMN_KO_vs_WT)
results_edgeR_D16_BMN_KO_vs_WT <- topTags(et_D16_BMN_KO_vs_WT, n = nrow(x), sort.by = "none")
head(results_edgeR_D16_BMN_KO_vs_WT$table)
edgeRt_D16_BMN_KO_vs_WT <- results_edgeR_D16_BMN_KO_vs_WT$table

z_D16_DAC_BMN_KO_vs_WT <- estimateDisp(D16_DAC_BMN_KO_vs_WT)
print(paste("sqrt of common dispersion (BCV) of","D16_DAC_BMN_KO", "Vs", "D16_DAC_BMN_WT", "is", sqrt(z_D16_DAC_BMN_KO_vs_WT$common.dispersion)))
#[1] "sqrt of common dispersion (BCV) of D16_DAC_BMN_KO Vs D16_DAC_BMN_WT is 0.314547116169699"
et_D16_DAC_BMN_KO_vs_WT <- exactTest(z_D16_DAC_BMN_KO_vs_WT)
results_edgeR_D16_DAC_BMN_KO_vs_WT <- topTags(et_D16_DAC_BMN_KO_vs_WT, n = nrow(x), sort.by = "none")
head(results_edgeR_D16_DAC_BMN_KO_vs_WT$table)
edgeRt_D16_DAC_BMN_KO_vs_WT <- results_edgeR_D16_DAC_BMN_KO_vs_WT$table

z_D26_DAC_KO_vs_WT <- estimateDisp(D26_DAC_KO_vs_WT)
print(paste("sqrt of common dispersion (BCV) of","D26_DAC_KO", "Vs", "D26_DAC_WT", "is", sqrt(z_D26_DAC_KO_vs_WT$common.dispersion)))
#[1] "sqrt of common dispersion (BCV) of D26_DAC_KO Vs D26_DAC_WT is 0.332243461458857"
et_D26_DAC_KO_vs_WT <- exactTest(z_D26_DAC_KO_vs_WT)
results_edgeR_D26_DAC_KO_vs_WT <- topTags(et_D26_DAC_KO_vs_WT, n = nrow(x), sort.by = "none")
head(results_edgeR_D26_DAC_KO_vs_WT$table)
edgeRt_D26_DAC_KO_vs_WT <- results_edgeR_D26_DAC_KO_vs_WT$table

z_D26_BMN_KO_vs_WT <- estimateDisp(D26_BMN_KO_vs_WT)
print(paste("sqrt of common dispersion (BCV) of","D26_BMN_KO", "Vs", "D26_BMN_WT", "is", sqrt(z_D26_BMN_KO_vs_WT$common.dispersion)))
#[1] "sqrt of common dispersion (BCV) of D26_BMN_KO Vs D26_BMN_WT is 0.346683211568326"
et_D26_BMN_KO_vs_WT <- exactTest(z_D26_BMN_KO_vs_WT)
results_edgeR_D26_BMN_KO_vs_WT <- topTags(et_D26_BMN_KO_vs_WT, n = nrow(x), sort.by = "none")
head(results_edgeR_D26_BMN_KO_vs_WT$table)
edgeRt_D26_BMN_KO_vs_WT <- results_edgeR_D26_BMN_KO_vs_WT$table

z_D26_DAC_BMN_KO_vs_WT <- estimateDisp(D26_DAC_BMN_KO_vs_WT)
print(paste("sqrt of common dispersion (BCV) of","D26_DAC_BMN_KO", "Vs", "D26_DAC_BMN_WT", "is", sqrt(z_D26_DAC_BMN_KO_vs_WT$common.dispersion)))
#[1] "sqrt of common dispersion (BCV) of D26_DAC_BMN_KO Vs D26_DAC_BMN_WT is 0.326485132870974"
et_D26_DAC_BMN_KO_vs_WT <- exactTest(z_D26_DAC_BMN_KO_vs_WT)
results_edgeR_D26_DAC_BMN_KO_vs_WT <- topTags(et_D26_DAC_BMN_KO_vs_WT, n = nrow(x), sort.by = "none")
head(results_edgeR_D26_DAC_BMN_KO_vs_WT$table)
edgeRt_D26_DAC_BMN_KO_vs_WT <- results_edgeR_D26_DAC_BMN_KO_vs_WT$table

#FDR<0.05  |logFC| >1
filter_FDR_0.05_NT_KO_vs_WT_up <- edgeRt_NT_KO_vs_WT[edgeRt_NT_KO_vs_WT$FDR < 0.05 & edgeRt_NT_KO_vs_WT$logFC > 1,]
filter_FDR_0.05_NT_KO_vs_WT_down <- edgeRt_NT_KO_vs_WT[edgeRt_NT_KO_vs_WT$FDR < 0.05 & edgeRt_NT_KO_vs_WT$logFC < -1,]
filter_FDR_0.05_D5_DAC_KO_vs_WT_up <- edgeRt_D5_DAC_KO_vs_WT[edgeRt_D5_DAC_KO_vs_WT$FDR < 0.05 & edgeRt_D5_DAC_KO_vs_WT$logFC > 1,]
filter_FDR_0.05_D5_DAC_KO_vs_WT_down <- edgeRt_D5_DAC_KO_vs_WT[edgeRt_D5_DAC_KO_vs_WT$FDR < 0.05 & edgeRt_D5_DAC_KO_vs_WT$logFC < -1,]
filter_FDR_0.05_D5_BMN_KO_vs_WT_up <- edgeRt_D5_BMN_KO_vs_WT[edgeRt_D5_BMN_KO_vs_WT$FDR < 0.05 & edgeRt_D5_BMN_KO_vs_WT$logFC > 1,]
filter_FDR_0.05_D5_BMN_KO_vs_WT_down <- edgeRt_D5_BMN_KO_vs_WT[edgeRt_D5_BMN_KO_vs_WT$FDR < 0.05 & edgeRt_D5_BMN_KO_vs_WT$logFC < -1,]
filter_FDR_0.05_D5_DAC_BMN_KO_vs_WT_up <- edgeRt_D5_DAC_BMN_KO_vs_WT[edgeRt_D5_DAC_BMN_KO_vs_WT$FDR < 0.05 & edgeRt_D5_DAC_BMN_KO_vs_WT$logFC > 1,]
filter_FDR_0.05_D5_DAC_BMN_KO_vs_WT_down <- edgeRt_D5_DAC_BMN_KO_vs_WT[edgeRt_D5_DAC_BMN_KO_vs_WT$FDR < 0.05 & edgeRt_D5_DAC_BMN_KO_vs_WT$logFC < -1,]
filter_FDR_0.05_NT_KO_vs_WT_up <- edgeRt_NT_KO_vs_WT[edgeRt_NT_KO_vs_WT$FDR < 0.05 & edgeRt_NT_KO_vs_WT$logFC > 1,]
filter_FDR_0.05_NT_KO_vs_WT_down <- edgeRt_NT_KO_vs_WT[edgeRt_NT_KO_vs_WT$FDR < 0.05 & edgeRt_NT_KO_vs_WT$logFC < -1,]
filter_FDR_0.05_D16_DAC_KO_vs_WT_up <- edgeRt_D16_DAC_KO_vs_WT[edgeRt_D16_DAC_KO_vs_WT$FDR < 0.05 & edgeRt_D16_DAC_KO_vs_WT$logFC > 1,]
filter_FDR_0.05_D16_DAC_KO_vs_WT_down <- edgeRt_D16_DAC_KO_vs_WT[edgeRt_D16_DAC_KO_vs_WT$FDR < 0.05 & edgeRt_D16_DAC_KO_vs_WT$logFC < -1,]
filter_FDR_0.05_D16_BMN_KO_vs_WT_up <- edgeRt_D16_BMN_KO_vs_WT[edgeRt_D16_BMN_KO_vs_WT$FDR < 0.05 & edgeRt_D16_BMN_KO_vs_WT$logFC > 1,]
filter_FDR_0.05_D16_BMN_KO_vs_WT_down <- edgeRt_D16_BMN_KO_vs_WT[edgeRt_D16_BMN_KO_vs_WT$FDR < 0.05 & edgeRt_D16_BMN_KO_vs_WT$logFC < -1,]
filter_FDR_0.05_D16_DAC_BMN_KO_vs_WT_up <- edgeRt_D16_DAC_BMN_KO_vs_WT[edgeRt_D16_DAC_BMN_KO_vs_WT$FDR < 0.05 & edgeRt_D16_DAC_BMN_KO_vs_WT$logFC > 1,]
filter_FDR_0.05_D16_DAC_BMN_KO_vs_WT_down <- edgeRt_D16_DAC_BMN_KO_vs_WT[edgeRt_D16_DAC_BMN_KO_vs_WT$FDR < 0.05 & edgeRt_D16_DAC_BMN_KO_vs_WT$logFC < -1,]
filter_FDR_0.05_NT_KO_vs_WT_up <- edgeRt_NT_KO_vs_WT[edgeRt_NT_KO_vs_WT$FDR < 0.05 & edgeRt_NT_KO_vs_WT$logFC > 1,]
filter_FDR_0.05_NT_KO_vs_WT_down <- edgeRt_NT_KO_vs_WT[edgeRt_NT_KO_vs_WT$FDR < 0.05 & edgeRt_NT_KO_vs_WT$logFC < -1,]
filter_FDR_0.05_D26_DAC_KO_vs_WT_up <- edgeRt_D26_DAC_KO_vs_WT[edgeRt_D26_DAC_KO_vs_WT$FDR < 0.05 & edgeRt_D26_DAC_KO_vs_WT$logFC > 1,]
filter_FDR_0.05_D26_DAC_KO_vs_WT_down <- edgeRt_D26_DAC_KO_vs_WT[edgeRt_D26_DAC_KO_vs_WT$FDR < 0.05 & edgeRt_D26_DAC_KO_vs_WT$logFC < -1,]
filter_FDR_0.05_D26_BMN_KO_vs_WT_up <- edgeRt_D26_BMN_KO_vs_WT[edgeRt_D26_BMN_KO_vs_WT$FDR < 0.05 & edgeRt_D26_BMN_KO_vs_WT$logFC > 1,]
filter_FDR_0.05_D26_BMN_KO_vs_WT_down <- edgeRt_D26_BMN_KO_vs_WT[edgeRt_D26_BMN_KO_vs_WT$FDR < 0.05 & edgeRt_D26_BMN_KO_vs_WT$logFC < -1,]
filter_FDR_0.05_D26_DAC_BMN_KO_vs_WT_up <- edgeRt_D26_DAC_BMN_KO_vs_WT[edgeRt_D26_DAC_BMN_KO_vs_WT$FDR < 0.05 & edgeRt_D26_DAC_BMN_KO_vs_WT$logFC > 1,]
filter_FDR_0.05_D26_DAC_BMN_KO_vs_WT_down <- edgeRt_D26_DAC_BMN_KO_vs_WT[edgeRt_D26_DAC_BMN_KO_vs_WT$FDR < 0.05 & edgeRt_D26_DAC_BMN_KO_vs_WT$logFC < -1,]

#save list
write.table(SINE,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/SINE.txt",row.names = T,col.names = T,quote = F)
write.table(LINE,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/LINE.txt",row.names = T,col.names = T,quote = F)
write.table(LTR,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/LTR.txt",row.names = T,col.names = T,quote = F)
write.table(z,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/z.txt",row.names = T,col.names = T,quote = F)
#UP
write.table(filter_FDR_0.05_NT_KO_vs_WT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/NT_KO_vs_WT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_D5_DAC_KO_vs_WT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D5_DAC_KO_vs_WT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_D5_BMN_KO_vs_WT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D5_BMN_KO_vs_WT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_D5_DAC_BMN_KO_vs_WT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D5_DAC_BMN_KO_vs_WT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_D16_DAC_KO_vs_WT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D16_DAC_KO_vs_WT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_D16_BMN_KO_vs_WT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D16_BMN_KO_vs_WT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_D16_DAC_BMN_KO_vs_WT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D16_DAC_BMN_KO_vs_WT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_D26_DAC_KO_vs_WT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D26_DAC_KO_vs_WT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_D26_BMN_KO_vs_WT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D26_BMN_KO_vs_WT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_D26_DAC_BMN_KO_vs_WT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D26_DAC_BMN_KO_vs_WT_up.txt",row.names = T,col.names = T,quote = F)

#DOWN
write.table(filter_FDR_0.05_NT_KO_vs_WT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/NT_KO_vs_WT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_D5_DAC_KO_vs_WT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D5_DAC_KO_vs_WT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_D5_BMN_KO_vs_WT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D5_BMN_KO_vs_WT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_D5_DAC_BMN_KO_vs_WT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D5_DAC_BMN_KO_vs_WT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_D16_DAC_KO_vs_WT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D16_DAC_KO_vs_WT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_D16_BMN_KO_vs_WT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D16_BMN_KO_vs_WT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_D16_DAC_BMN_KO_vs_WT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D16_DAC_BMN_KO_vs_WT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_D26_DAC_KO_vs_WT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D26_DAC_KO_vs_WT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_D26_BMN_KO_vs_WT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D26_BMN_KO_vs_WT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_D26_DAC_BMN_KO_vs_WT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/D26_DAC_BMN_KO_vs_WT_down.txt",row.names = T,col.names = T,quote = F)


#KO vs WT
#TE
TE_up_KO_vs_WT_DAC <- c(row.names(D5_DAC_KO_vs_WT_up),row.names(D16_DAC_KO_vs_WT_up),row.names(D26_DAC_KO_vs_WT_up));TE_up_KO_vs_WT_DAC <- transform(TE_up_KO_vs_WT_DAC,duplication=duplicated(TE_up_KO_vs_WT_DAC))
TE_up_KO_vs_WT_DAC <- TE_up_KO_vs_WT_DAC[TE_up_KO_vs_WT_DAC$duplication==F,];row.names(TE_up_KO_vs_WT_DAC) <- TE_up_KO_vs_WT_DAC$X_data
TE_up_KO_vs_WT_DAC_excludeNT <- transform(TE_up_KO_vs_WT_DAC,isinNT=row.names(TE_up_KO_vs_WT_DAC)%in%row.names(NT_KO_vs_WT_up));TE_up_KO_vs_WT_DAC_excludeNT <- TE_up_KO_vs_WT_DAC_excludeNT[TE_up_KO_vs_WT_DAC_excludeNT$isinNT==F,]

TE_up_KO_vs_WT_BMN <- c(row.names(D5_BMN_KO_vs_WT_up),row.names(D16_BMN_KO_vs_WT_up),row.names(D26_BMN_KO_vs_WT_up));TE_up_KO_vs_WT_BMN <- transform(TE_up_KO_vs_WT_BMN,duplication=duplicated(TE_up_KO_vs_WT_BMN))
TE_up_KO_vs_WT_BMN <- TE_up_KO_vs_WT_BMN[TE_up_KO_vs_WT_BMN$duplication==F,];row.names(TE_up_KO_vs_WT_BMN) <- TE_up_KO_vs_WT_BMN$X_data
TE_up_KO_vs_WT_BMN_excludeNT <- transform(TE_up_KO_vs_WT_BMN,isinNT=row.names(TE_up_KO_vs_WT_BMN)%in%row.names(NT_KO_vs_WT_up));TE_up_KO_vs_WT_BMN_excludeNT <- TE_up_KO_vs_WT_BMN_excludeNT[TE_up_KO_vs_WT_BMN_excludeNT$isinNT==F,]

TE_up_KO_vs_WT_DAC_BMN <- c(row.names(D5_DAC_BMN_KO_vs_WT_up),row.names(D16_DAC_BMN_KO_vs_WT_up),row.names(D26_DAC_BMN_KO_vs_WT_up));TE_up_KO_vs_WT_DAC_BMN <- transform(TE_up_KO_vs_WT_DAC_BMN,duplication=duplicated(TE_up_KO_vs_WT_DAC_BMN))
TE_up_KO_vs_WT_DAC_BMN <- TE_up_KO_vs_WT_DAC_BMN[TE_up_KO_vs_WT_DAC_BMN$duplication==F,];row.names(TE_up_KO_vs_WT_DAC_BMN) <- TE_up_KO_vs_WT_DAC_BMN$X_data
TE_up_KO_vs_WT_DAC_BMN_excludeNT <- transform(TE_up_KO_vs_WT_DAC_BMN,isinNT=row.names(TE_up_KO_vs_WT_DAC_BMN)%in%row.names(NT_KO_vs_WT_up));TE_up_KO_vs_WT_DAC_BMN_excludeNT <- TE_up_KO_vs_WT_DAC_BMN_excludeNT[TE_up_KO_vs_WT_DAC_BMN_excludeNT$isinNT==F,]

TE_up_KO_vs_WT_ALL <- c(row.names(D5_DAC_KO_vs_WT_up),row.names(D16_DAC_KO_vs_WT_up),row.names(D26_DAC_KO_vs_WT_up),row.names(D5_BMN_KO_vs_WT_up),row.names(D16_BMN_KO_vs_WT_up),row.names(D26_BMN_KO_vs_WT_up),row.names(D5_DAC_BMN_KO_vs_WT_up),row.names(D16_DAC_BMN_KO_vs_WT_up),row.names(D26_DAC_BMN_KO_vs_WT_up));TE_up_KO_vs_WT_ALL <- transform(TE_up_KO_vs_WT_ALL,duplication=duplicated(TE_up_KO_vs_WT_ALL))
TE_up_KO_vs_WT_ALL <- TE_up_KO_vs_WT_ALL[TE_up_KO_vs_WT_ALL$duplication==F,];row.names(TE_up_KO_vs_WT_ALL) <- TE_up_KO_vs_WT_ALL$X_data
TE_up_KO_vs_WT_ALL_excludeNT <- transform(TE_up_KO_vs_WT_ALL,isinNT=row.names(TE_up_KO_vs_WT_ALL)%in%row.names(NT_KO_vs_WT_up));TE_up_KO_vs_WT_ALL_excludeNT <- TE_up_KO_vs_WT_ALL_excludeNT[TE_up_KO_vs_WT_ALL_excludeNT$isinNT==F,]

