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

###HEATMAP
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

###counts
combined_x <- matrix(nr=551005,nc=20)
combined_x[,1] <- rowMeans(x[,c(1,2)])
combined_x[,2] <- rowMeans(x[,c(3,4)])
combined_x[,3] <- rowMeans(x[,c(5,6)])
combined_x[,4] <- rowMeans(x[,c(7,8)])
combined_x[,5] <- rowMeans(x[,c(9,10)])
combined_x[,6] <- rowMeans(x[,c(11,12)])
combined_x[,7] <- rowMeans(x[,c(13,14)])
combined_x[,8] <- rowMeans(x[,c(15,16)])
combined_x[,9] <- rowMeans(x[,c(17,18)])
combined_x[,10] <- rowMeans(x[,c(19,20)])
combined_x[,11] <- rowMeans(x[,c(21,22)])
combined_x[,12] <- rowMeans(x[,c(23,24)])
combined_x[,13] <- rowMeans(x[,c(25,26)])
combined_x[,14] <- rowMeans(x[,c(27,28)])
combined_x[,15] <- rowMeans(x[,c(29,30)])
combined_x[,16] <- rowMeans(x[,c(31,32)])
combined_x[,17] <- rowMeans(x[,c(33,34)])
combined_x[,18] <- rowMeans(x[,c(35,36)])
combined_x[,19] <- rowMeans(x[,c(37,38)])
combined_x[,20] <- rowMeans(x[,c(39,40)])
colnames(combined_x) <- c("WT_NT","Day5_WT_BMN","Day16_WT_BMN","Day26_WT_BMN","Day5_WT_DAC","Day16_WT_DAC","Day26_WT_DAC","Day5_WT_DAC+BMN","Day16_WT_DAC+BMN","Day26_WT_DAC+BMN","KO_NT","Day5_KO_BMN","Day16_KO_BMN","Day26_KO_BMN","Day5_KO_DAC","Day16_KO_DAC","Day26_KO_DAC","Day5_KO_DAC+BMN","Day16_KO_DAC+BMN","Day26_KO_DAC+BMN")
rownames(combined_x) <- rownames(x)

TEsum <- matrix(colSums(combined_x))
rownames(TEsum) <- colnames(combined_x)
write.table(TEsum,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/TEsum.txt",row.names = T,col.names = T,quote = F,sep = "\t")
