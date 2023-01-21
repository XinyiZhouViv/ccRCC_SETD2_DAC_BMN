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
#先算filter基因后的z-score
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
#未经筛选的SINE, LINE, LTR list
SINE <- myTEdata[myTEdata$Classid=="SINE",c(1,2)]
LINE <- myTEdata[myTEdata$Classid=="LINE",c(1,2)]
LTR <- myTEdata[myTEdata$Classid=="LTR",c(1,2)]

#将y分成KO vs WT的各组进行比较
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

#将y分成Treat vs NT的各组进行比较
WT_D5_DAC_vs_NT <- y[,c(1,2,9,10)]
WT_D5_BMN_vs_NT <- y[,c(1,2,3,4)]
WT_D5_DAC_BMN_vs_NT <- y[,c(1,2,15,16)]
WT_D16_DAC_vs_NT <- y[,c(1,2,11,12)]
WT_D16_BMN_vs_NT <- y[,c(1,2,5,6)]
WT_D16_DAC_BMN_vs_NT <- y[,c(1,2,17,18)]
WT_D26_DAC_vs_NT <- y[,c(1,2,13,14)]
WT_D26_BMN_vs_NT <- y[,c(1,2,7,8)]
WT_D26_DAC_BMN_vs_NT <- y[,c(1,2,19,20)]
KO_D5_DAC_vs_NT <- y[,c(21,22,29,30)]
KO_D5_BMN_vs_NT <- y[,c(21,22,23,24)]
KO_D5_DAC_BMN_vs_NT <- y[,c(21,22,35,36)]
KO_D16_DAC_vs_NT <- y[,c(21,22,31,32)]
KO_D16_BMN_vs_NT <- y[,c(21,22,25,26)]
KO_D16_DAC_BMN_vs_NT <-y[,c(21,22,37,38)]
KO_D26_DAC_vs_NT <- y[,c(21,22,33,34)]
KO_D26_BMN_vs_NT <- y[,c(21,22,27,28)]
KO_D26_DAC_BMN_vs_NT <- y[,c(21,22,39,40)]

z_WT_D5_DAC_vs_NT <- estimateDisp(WT_D5_DAC_vs_NT)
print(paste("sqrt of common dispersion (BCV) of WT Day5","DAC", "Vs", "NT", "is", sqrt(z_WT_D5_DAC_vs_NT$common.dispersion)))
# [1] "sqrt of common dispersion (BCV) of WT Day5 DAC Vs NT is 0.368956813227943"
et_WT_D5_DAC_vs_NT <- exactTest(z_WT_D5_DAC_vs_NT)
results_edgeR_WT_D5_DAC_vs_NT <- topTags(et_WT_D5_DAC_vs_NT, n = nrow(x), sort.by = "none")
head(results_edgeR_WT_D5_DAC_vs_NT$table)
edgeRt_WT_D5_DAC_vs_NT <- results_edgeR_WT_D5_DAC_vs_NT$table

z_WT_D5_BMN_vs_NT <- estimateDisp(WT_D5_BMN_vs_NT)
print(paste("sqrt of common dispersion (BCV) of WT Day5","BMN", "Vs", "NT", "is", sqrt(z_WT_D5_BMN_vs_NT$common.dispersion)))
# [1] "sqrt of common dispersion (BCV) of WT Day5 BMN Vs NT is 0.357096861043914"
et_WT_D5_BMN_vs_NT <- exactTest(z_WT_D5_BMN_vs_NT)
results_edgeR_WT_D5_BMN_vs_NT <- topTags(et_WT_D5_BMN_vs_NT, n = nrow(x), sort.by = "none")
head(results_edgeR_WT_D5_BMN_vs_NT$table)
edgeRt_WT_D5_BMN_vs_NT <- results_edgeR_WT_D5_BMN_vs_NT$table

z_WT_D5_DAC_BMN_vs_NT <- estimateDisp(WT_D5_DAC_BMN_vs_NT)
print(paste("sqrt of common dispersion (BCV) of WT Day5","DAC+BMN", "Vs", "NT", "is", sqrt(z_WT_D5_DAC_BMN_vs_NT$common.dispersion)))
# [1] "sqrt of common dispersion (BCV) of WT Day5 DAC+BMN Vs NT is 0.362846536762564"
et_WT_D5_DAC_BMN_vs_NT <- exactTest(z_WT_D5_DAC_BMN_vs_NT)
results_edgeR_WT_D5_DAC_BMN_vs_NT <- topTags(et_WT_D5_DAC_BMN_vs_NT, n = nrow(x), sort.by = "none")
head(results_edgeR_WT_D5_DAC_BMN_vs_NT$table)
edgeRt_WT_D5_DAC_BMN_vs_NT <- results_edgeR_WT_D5_DAC_BMN_vs_NT$table

z_WT_D16_DAC_vs_NT <- estimateDisp(WT_D16_DAC_vs_NT)
print(paste("sqrt of common dispersion (BCV) of WT Day16","DAC", "Vs", "NT", "is", sqrt(z_WT_D16_DAC_vs_NT$common.dispersion)))
# [1] "sqrt of common dispersion (BCV) of WT Day16 DAC Vs NT is 0.34653776288664"
et_WT_D16_DAC_vs_NT <- exactTest(z_WT_D16_DAC_vs_NT)
results_edgeR_WT_D16_DAC_vs_NT <- topTags(et_WT_D16_DAC_vs_NT, n = nrow(x), sort.by = "none")
head(results_edgeR_WT_D16_DAC_vs_NT$table)
edgeRt_WT_D16_DAC_vs_NT <- results_edgeR_WT_D16_DAC_vs_NT$table

z_WT_D16_BMN_vs_NT <- estimateDisp(WT_D16_BMN_vs_NT)
print(paste("sqrt of common dispersion (BCV) of WT Day16","BMN", "Vs", "NT", "is", sqrt(z_WT_D16_BMN_vs_NT$common.dispersion)))
# [1] "sqrt of common dispersion (BCV) of WT Day16 BMN Vs NT is 0.357504956879879"
et_WT_D16_BMN_vs_NT <- exactTest(z_WT_D16_BMN_vs_NT)
results_edgeR_WT_D16_BMN_vs_NT <- topTags(et_WT_D16_BMN_vs_NT, n = nrow(x), sort.by = "none")
head(results_edgeR_WT_D16_BMN_vs_NT$table)
edgeRt_WT_D16_BMN_vs_NT <- results_edgeR_WT_D16_BMN_vs_NT$table

z_WT_D16_DAC_BMN_vs_NT <- estimateDisp(WT_D16_DAC_BMN_vs_NT)
print(paste("sqrt of common dispersion (BCV) of WT Day16","DAC+BMN", "Vs", "NT", "is", sqrt(z_WT_D16_DAC_BMN_vs_NT$common.dispersion)))
# [1] "sqrt of common dispersion (BCV) of WT Day16 DAC+BMN Vs NT is 0.354638705557314"
et_WT_D16_DAC_BMN_vs_NT <- exactTest(z_WT_D16_DAC_BMN_vs_NT)
results_edgeR_WT_D16_DAC_BMN_vs_NT <- topTags(et_WT_D16_DAC_BMN_vs_NT, n = nrow(x), sort.by = "none")
head(results_edgeR_WT_D16_DAC_BMN_vs_NT$table)
edgeRt_WT_D16_DAC_BMN_vs_NT <- results_edgeR_WT_D16_DAC_BMN_vs_NT$table

z_WT_D26_DAC_vs_NT <- estimateDisp(WT_D26_DAC_vs_NT)
print(paste("sqrt of common dispersion (BCV) of WT Day26","DAC", "Vs", "NT", "is", sqrt(z_WT_D26_DAC_vs_NT$common.dispersion)))
# [1] "sqrt of common dispersion (BCV) of WT Day26 DAC Vs NT is 0.390073276623852"
et_WT_D26_DAC_vs_NT <- exactTest(z_WT_D26_DAC_vs_NT)
results_edgeR_WT_D26_DAC_vs_NT <- topTags(et_WT_D26_DAC_vs_NT, n = nrow(x), sort.by = "none")
head(results_edgeR_WT_D26_DAC_vs_NT$table)
edgeRt_WT_D26_DAC_vs_NT <- results_edgeR_WT_D26_DAC_vs_NT$table

z_WT_D26_BMN_vs_NT <- estimateDisp(WT_D26_BMN_vs_NT)
print(paste("sqrt of common dispersion (BCV) of WT Day26","BMN", "Vs", "NT", "is", sqrt(z_WT_D26_BMN_vs_NT$common.dispersion)))
# [1] "sqrt of common dispersion (BCV) of WT Day26 BMN Vs NT is 0.377530274631764"
et_WT_D26_BMN_vs_NT <- exactTest(z_WT_D26_BMN_vs_NT)
results_edgeR_WT_D26_BMN_vs_NT <- topTags(et_WT_D26_BMN_vs_NT, n = nrow(x), sort.by = "none")
head(results_edgeR_WT_D26_BMN_vs_NT$table)
edgeRt_WT_D26_BMN_vs_NT <- results_edgeR_WT_D26_BMN_vs_NT$table

z_WT_D26_DAC_BMN_vs_NT <- estimateDisp(WT_D26_DAC_BMN_vs_NT)
print(paste("sqrt of common dispersion (BCV) of WT Day26","DAC+BMN", "Vs", "NT", "is", sqrt(z_WT_D26_DAC_BMN_vs_NT$common.dispersion)))
# [1] "sqrt of common dispersion (BCV) of WT Day26 DAC+BMN Vs NT is 0.365728638477629"
et_WT_D26_DAC_BMN_vs_NT <- exactTest(z_WT_D26_DAC_BMN_vs_NT)
results_edgeR_WT_D26_DAC_BMN_vs_NT <- topTags(et_WT_D26_DAC_BMN_vs_NT, n = nrow(x), sort.by = "none")
head(results_edgeR_WT_D26_DAC_BMN_vs_NT$table)
edgeRt_WT_D26_DAC_BMN_vs_NT <- results_edgeR_WT_D26_DAC_BMN_vs_NT$table

z_KO_D5_DAC_vs_NT <- estimateDisp(KO_D5_DAC_vs_NT)
print(paste("sqrt of common dispersion (BCV) of KO Day5","DAC", "Vs", "NT", "is", sqrt(z_KO_D5_DAC_vs_NT$common.dispersion)))
# [1] "sqrt of common dispersion (BCV) of KO Day5 DAC Vs NT is 0.353748772788042"
et_KO_D5_DAC_vs_NT <- exactTest(z_KO_D5_DAC_vs_NT)
results_edgeR_KO_D5_DAC_vs_NT <- topTags(et_KO_D5_DAC_vs_NT, n = nrow(x), sort.by = "none")
head(results_edgeR_KO_D5_DAC_vs_NT$table)
edgeRt_KO_D5_DAC_vs_NT <- results_edgeR_KO_D5_DAC_vs_NT$table

z_KO_D5_BMN_vs_NT <- estimateDisp(KO_D5_BMN_vs_NT)
print(paste("sqrt of common dispersion (BCV) of KO Day5","BMN", "Vs", "NT", "is", sqrt(z_KO_D5_BMN_vs_NT$common.dispersion)))
# [1] "sqrt of common dispersion (BCV) of KO Day5 BMN Vs NT is 0.338879963209856"
et_KO_D5_BMN_vs_NT <- exactTest(z_KO_D5_BMN_vs_NT)
results_edgeR_KO_D5_BMN_vs_NT <- topTags(et_KO_D5_BMN_vs_NT, n = nrow(x), sort.by = "none")
head(results_edgeR_KO_D5_BMN_vs_NT$table)
edgeRt_KO_D5_BMN_vs_NT <- results_edgeR_KO_D5_BMN_vs_NT$table

z_KO_D5_DAC_BMN_vs_NT <- estimateDisp(KO_D5_DAC_BMN_vs_NT)
print(paste("sqrt of common dispersion (BCV) of KO Day5","DAC+BMN", "Vs", "NT", "is", sqrt(z_KO_D5_DAC_BMN_vs_NT$common.dispersion)))
# [1] "sqrt of common dispersion (BCV) of KO Day5 DAC+BMN Vs NT is 0.314944636355001"
et_KO_D5_DAC_BMN_vs_NT <- exactTest(z_KO_D5_DAC_BMN_vs_NT)
results_edgeR_KO_D5_DAC_BMN_vs_NT <- topTags(et_KO_D5_DAC_BMN_vs_NT, n = nrow(x), sort.by = "none")
head(results_edgeR_KO_D5_DAC_BMN_vs_NT$table)
edgeRt_KO_D5_DAC_BMN_vs_NT <- results_edgeR_KO_D5_DAC_BMN_vs_NT$table

z_KO_D16_DAC_vs_NT <- estimateDisp(KO_D16_DAC_vs_NT)
print(paste("sqrt of common dispersion (BCV) of KO Day16","DAC", "Vs", "NT", "is", sqrt(z_KO_D16_DAC_vs_NT$common.dispersion)))
# [1] "sqrt of common dispersion (BCV) of KO Day16 DAC Vs NT is 0.330929857059437"
et_KO_D16_DAC_vs_NT <- exactTest(z_KO_D16_DAC_vs_NT)
results_edgeR_KO_D16_DAC_vs_NT <- topTags(et_KO_D16_DAC_vs_NT, n = nrow(x), sort.by = "none")
head(results_edgeR_KO_D16_DAC_vs_NT$table)
edgeRt_KO_D16_DAC_vs_NT <- results_edgeR_KO_D16_DAC_vs_NT$table

z_KO_D16_BMN_vs_NT <- estimateDisp(KO_D16_BMN_vs_NT)
print(paste("sqrt of common dispersion (BCV) of KO Day16","BMN", "Vs", "NT", "is", sqrt(z_KO_D16_BMN_vs_NT$common.dispersion)))
# [1] "sqrt of common dispersion (BCV) of KO Day16 BMN Vs NT is 0.326931720300696"
et_KO_D16_BMN_vs_NT <- exactTest(z_KO_D16_BMN_vs_NT)
results_edgeR_KO_D16_BMN_vs_NT <- topTags(et_KO_D16_BMN_vs_NT, n = nrow(x), sort.by = "none")
head(results_edgeR_KO_D16_BMN_vs_NT$table)
edgeRt_KO_D16_BMN_vs_NT <- results_edgeR_KO_D16_BMN_vs_NT$table

z_KO_D16_DAC_BMN_vs_NT <- estimateDisp(KO_D16_DAC_BMN_vs_NT)
print(paste("sqrt of common dispersion (BCV) of KO Day16","DAC+BMN", "Vs", "NT", "is", sqrt(z_KO_D16_DAC_BMN_vs_NT$common.dispersion)))
# [1] "sqrt of common dispersion (BCV) of KO Day16 DAC+BMN Vs NT is 0.326395743579499"
et_KO_D16_DAC_BMN_vs_NT <- exactTest(z_KO_D16_DAC_BMN_vs_NT)
results_edgeR_KO_D16_DAC_BMN_vs_NT <- topTags(et_KO_D16_DAC_BMN_vs_NT, n = nrow(x), sort.by = "none")
head(results_edgeR_KO_D16_DAC_BMN_vs_NT$table)
edgeRt_KO_D16_DAC_BMN_vs_NT <- results_edgeR_KO_D16_DAC_BMN_vs_NT$table

z_KO_D26_DAC_vs_NT <- estimateDisp(KO_D26_DAC_vs_NT)
print(paste("sqrt of common dispersion (BCV) of KO Day26","DAC", "Vs", "NT", "is", sqrt(z_KO_D26_DAC_vs_NT$common.dispersion)))
# "sqrt of common dispersion (BCV) of KO Day26 DAC Vs NT is 0.311368775439774"
et_KO_D26_DAC_vs_NT <- exactTest(z_KO_D26_DAC_vs_NT)
results_edgeR_KO_D26_DAC_vs_NT <- topTags(et_KO_D26_DAC_vs_NT, n = nrow(x), sort.by = "none")
head(results_edgeR_KO_D26_DAC_vs_NT$table)
edgeRt_KO_D26_DAC_vs_NT <- results_edgeR_KO_D26_DAC_vs_NT$table

z_KO_D26_BMN_vs_NT <- estimateDisp(KO_D26_BMN_vs_NT)
print(paste("sqrt of common dispersion (BCV) of KO Day26","BMN", "Vs", "NT", "is", sqrt(z_KO_D26_BMN_vs_NT$common.dispersion)))
# [1] "sqrt of common dispersion (BCV) of KO Day26 BMN Vs NT is 0.33371195520111"
et_KO_D26_BMN_vs_NT <- exactTest(z_KO_D26_BMN_vs_NT)
results_edgeR_KO_D26_BMN_vs_NT <- topTags(et_KO_D26_BMN_vs_NT, n = nrow(x), sort.by = "none")
head(results_edgeR_KO_D26_BMN_vs_NT$table)
edgeRt_KO_D26_BMN_vs_NT <- results_edgeR_KO_D26_BMN_vs_NT$table

z_KO_D26_DAC_BMN_vs_NT <- estimateDisp(KO_D26_DAC_BMN_vs_NT)
print(paste("sqrt of common dispersion (BCV) of KO Day26","DAC+BMN", "Vs", "NT", "is", sqrt(z_KO_D26_DAC_BMN_vs_NT$common.dispersion)))
# [1] "sqrt of common dispersion (BCV) of KO Day26 DAC+BMN Vs NT is 0.330581673737279"
et_KO_D26_DAC_BMN_vs_NT <- exactTest(z_KO_D26_DAC_BMN_vs_NT)
results_edgeR_KO_D26_DAC_BMN_vs_NT <- topTags(et_KO_D26_DAC_BMN_vs_NT, n = nrow(x), sort.by = "none")
head(results_edgeR_KO_D26_DAC_BMN_vs_NT$table)
edgeRt_KO_D26_DAC_BMN_vs_NT <- results_edgeR_KO_D26_DAC_BMN_vs_NT$table

#筛选FDR<0.05 和 |logFC| >1的基因
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

filter_FDR_0.05_WT_D5_DAC_vs_NT_up <- edgeRt_WT_D5_DAC_vs_NT[edgeRt_WT_D5_DAC_vs_NT$FDR < 0.05 & edgeRt_WT_D5_DAC_vs_NT$logFC > 1,]
filter_FDR_0.05_WT_D5_DAC_vs_NT_down <- edgeRt_WT_D5_DAC_vs_NT[edgeRt_WT_D5_DAC_vs_NT$FDR < 0.05 & edgeRt_WT_D5_DAC_vs_NT$logFC < -1,]
#filter_FDR_0.05_WT_D5_BMN_vs_NT_up <- edgeRt_WT_D5_BMN_vs_NT[edgeRt_WT_D5_BMN_vs_NT$FDR < 0.05 & edgeRt_WT_D5_BMN_vs_NT$logFC > 1,]无可用数据
#filter_FDR_0.05_WT_D5_BMN_vs_NT_down <- edgeRt_WT_D5_BMN_vs_NT[edgeRt_WT_D5_BMN_vs_NT$FDR < 0.05 & edgeRt_WT_D5_BMN_vs_NT$logFC < -1,]无可用数据
filter_FDR_0.05_WT_D5_DAC_BMN_vs_NT_up <- edgeRt_WT_D5_DAC_BMN_vs_NT[edgeRt_WT_D5_DAC_BMN_vs_NT$FDR < 0.05 & edgeRt_WT_D5_DAC_BMN_vs_NT$logFC > 1,]
filter_FDR_0.05_WT_D5_DAC_BMN_vs_NT_down <- edgeRt_WT_D5_DAC_BMN_vs_NT[edgeRt_WT_D5_DAC_BMN_vs_NT$FDR < 0.05 & edgeRt_WT_D5_DAC_BMN_vs_NT$logFC < -1,]
filter_FDR_0.05_WT_D16_DAC_vs_NT_up <- edgeRt_WT_D16_DAC_vs_NT[edgeRt_WT_D16_DAC_vs_NT$FDR < 0.05 & edgeRt_WT_D16_DAC_vs_NT$logFC > 1,]
filter_FDR_0.05_WT_D16_DAC_vs_NT_down <- edgeRt_WT_D16_DAC_vs_NT[edgeRt_WT_D16_DAC_vs_NT$FDR < 0.05 & edgeRt_WT_D16_DAC_vs_NT$logFC < -1,]
filter_FDR_0.05_WT_D16_BMN_vs_NT_up <- edgeRt_WT_D16_BMN_vs_NT[edgeRt_WT_D16_BMN_vs_NT$FDR < 0.05 & edgeRt_WT_D16_BMN_vs_NT$logFC > 1,]
#filter_FDR_0.05_WT_D16_BMN_vs_NT_down <- edgeRt_WT_D16_BMN_vs_NT[edgeRt_WT_D16_BMN_vs_NT$FDR < 0.05 & edgeRt_WT_D16_BMN_vs_NT$logFC < -1,]无可用数据
filter_FDR_0.05_WT_D16_DAC_BMN_vs_NT_up <- edgeRt_WT_D16_DAC_BMN_vs_NT[edgeRt_WT_D16_DAC_BMN_vs_NT$FDR < 0.05 & edgeRt_WT_D16_DAC_BMN_vs_NT$logFC > 1,]
filter_FDR_0.05_WT_D16_DAC_BMN_vs_NT_down <- edgeRt_WT_D16_DAC_BMN_vs_NT[edgeRt_WT_D16_DAC_BMN_vs_NT$FDR < 0.05 & edgeRt_WT_D16_DAC_BMN_vs_NT$logFC < -1,]
filter_FDR_0.05_WT_D26_DAC_vs_NT_up <- edgeRt_WT_D26_DAC_vs_NT[edgeRt_WT_D26_DAC_vs_NT$FDR < 0.05 & edgeRt_WT_D26_DAC_vs_NT$logFC > 1,]
filter_FDR_0.05_WT_D26_DAC_vs_NT_down <- edgeRt_WT_D26_DAC_vs_NT[edgeRt_WT_D26_DAC_vs_NT$FDR < 0.05 & edgeRt_WT_D26_DAC_vs_NT$logFC < -1,]
filter_FDR_0.05_WT_D26_BMN_vs_NT_up <- edgeRt_WT_D26_BMN_vs_NT[edgeRt_WT_D26_BMN_vs_NT$FDR < 0.05 & edgeRt_WT_D26_BMN_vs_NT$logFC > 1,]
filter_FDR_0.05_WT_D26_BMN_vs_NT_down <- edgeRt_WT_D26_BMN_vs_NT[edgeRt_WT_D26_BMN_vs_NT$FDR < 0.05 & edgeRt_WT_D26_BMN_vs_NT$logFC < -1,]
filter_FDR_0.05_WT_D26_DAC_BMN_vs_NT_up <- edgeRt_WT_D26_DAC_BMN_vs_NT[edgeRt_WT_D26_DAC_BMN_vs_NT$FDR < 0.05 & edgeRt_WT_D26_DAC_BMN_vs_NT$logFC > 1,]
filter_FDR_0.05_WT_D26_DAC_BMN_vs_NT_down <- edgeRt_WT_D26_DAC_BMN_vs_NT[edgeRt_WT_D26_DAC_BMN_vs_NT$FDR < 0.05 & edgeRt_WT_D26_DAC_BMN_vs_NT$logFC < -1,]

filter_FDR_0.05_KO_D5_DAC_vs_NT_up <- edgeRt_KO_D5_DAC_vs_NT[edgeRt_KO_D5_DAC_vs_NT$FDR < 0.05 & edgeRt_KO_D5_DAC_vs_NT$logFC > 1,]
filter_FDR_0.05_KO_D5_DAC_vs_NT_down <- edgeRt_KO_D5_DAC_vs_NT[edgeRt_KO_D5_DAC_vs_NT$FDR < 0.05 & edgeRt_KO_D5_DAC_vs_NT$logFC < -1,]
filter_FDR_0.05_KO_D5_BMN_vs_NT_up <- edgeRt_KO_D5_BMN_vs_NT[edgeRt_KO_D5_BMN_vs_NT$FDR < 0.05 & edgeRt_KO_D5_BMN_vs_NT$logFC > 1,]
filter_FDR_0.05_KO_D5_BMN_vs_NT_down <- edgeRt_KO_D5_BMN_vs_NT[edgeRt_KO_D5_BMN_vs_NT$FDR < 0.05 & edgeRt_KO_D5_BMN_vs_NT$logFC < -1,]
filter_FDR_0.05_KO_D5_DAC_BMN_vs_NT_up <- edgeRt_KO_D5_DAC_BMN_vs_NT[edgeRt_KO_D5_DAC_BMN_vs_NT$FDR < 0.05 & edgeRt_KO_D5_DAC_BMN_vs_NT$logFC > 1,]
filter_FDR_0.05_KO_D5_DAC_BMN_vs_NT_down <- edgeRt_KO_D5_DAC_BMN_vs_NT[edgeRt_KO_D5_DAC_BMN_vs_NT$FDR < 0.05 & edgeRt_KO_D5_DAC_BMN_vs_NT$logFC < -1,]
filter_FDR_0.05_KO_D16_DAC_vs_NT_up <- edgeRt_KO_D16_DAC_vs_NT[edgeRt_KO_D16_DAC_vs_NT$FDR < 0.05 & edgeRt_KO_D16_DAC_vs_NT$logFC > 1,]
filter_FDR_0.05_KO_D16_DAC_vs_NT_down <- edgeRt_KO_D16_DAC_vs_NT[edgeRt_KO_D16_DAC_vs_NT$FDR < 0.05 & edgeRt_KO_D16_DAC_vs_NT$logFC < -1,]
filter_FDR_0.05_KO_D16_BMN_vs_NT_up <- edgeRt_KO_D16_BMN_vs_NT[edgeRt_KO_D16_BMN_vs_NT$FDR < 0.05 & edgeRt_KO_D16_BMN_vs_NT$logFC > 1,]
filter_FDR_0.05_KO_D16_BMN_vs_NT_down <- edgeRt_KO_D16_BMN_vs_NT[edgeRt_KO_D16_BMN_vs_NT$FDR < 0.05 & edgeRt_KO_D16_BMN_vs_NT$logFC < -1,]
filter_FDR_0.05_KO_D16_DAC_BMN_vs_NT_up <- edgeRt_KO_D16_DAC_BMN_vs_NT[edgeRt_KO_D16_DAC_BMN_vs_NT$FDR < 0.05 & edgeRt_KO_D16_DAC_BMN_vs_NT$logFC > 1,]
filter_FDR_0.05_KO_D16_DAC_BMN_vs_NT_down <- edgeRt_KO_D16_DAC_BMN_vs_NT[edgeRt_KO_D16_DAC_BMN_vs_NT$FDR < 0.05 & edgeRt_KO_D16_DAC_BMN_vs_NT$logFC < -1,]
filter_FDR_0.05_KO_D26_DAC_vs_NT_up <- edgeRt_KO_D26_DAC_vs_NT[edgeRt_KO_D26_DAC_vs_NT$FDR < 0.05 & edgeRt_KO_D26_DAC_vs_NT$logFC > 1,]
filter_FDR_0.05_KO_D26_DAC_vs_NT_down <- edgeRt_KO_D26_DAC_vs_NT[edgeRt_KO_D26_DAC_vs_NT$FDR < 0.05 & edgeRt_KO_D26_DAC_vs_NT$logFC < -1,]
filter_FDR_0.05_KO_D26_BMN_vs_NT_up <- edgeRt_KO_D26_BMN_vs_NT[edgeRt_KO_D26_BMN_vs_NT$FDR < 0.05 & edgeRt_KO_D26_BMN_vs_NT$logFC > 1,]
filter_FDR_0.05_KO_D26_BMN_vs_NT_down <- edgeRt_KO_D26_BMN_vs_NT[edgeRt_KO_D26_BMN_vs_NT$FDR < 0.05 & edgeRt_KO_D26_BMN_vs_NT$logFC < -1,]
filter_FDR_0.05_KO_D26_DAC_BMN_vs_NT_up <- edgeRt_KO_D26_DAC_BMN_vs_NT[edgeRt_KO_D26_DAC_BMN_vs_NT$FDR < 0.05 & edgeRt_KO_D26_DAC_BMN_vs_NT$logFC > 1,]
filter_FDR_0.05_KO_D26_DAC_BMN_vs_NT_down <- edgeRt_KO_D26_DAC_BMN_vs_NT[edgeRt_KO_D26_DAC_BMN_vs_NT$FDR < 0.05 & edgeRt_KO_D26_DAC_BMN_vs_NT$logFC < -1,]

#保存总list，各组差异基因list
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
write.table(filter_FDR_0.05_WT_D5_DAC_vs_NT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/WT_D5_DAC_vs_NT_up.txt",row.names = T,col.names = T,quote = F)
#write.table(filter_FDR_0.05_WT_D5_BMN_vs_NT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/WT_D5_BMN_vs_NT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_WT_D5_DAC_BMN_vs_NT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/WT_D5_DAC_BMN_vs_NT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_WT_D16_DAC_vs_NT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/WT_D16_DAC_vs_NT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_WT_D16_BMN_vs_NT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/WT_D16_BMN_vs_NT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_WT_D16_DAC_BMN_vs_NT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/WT_D16_DAC_BMN_vs_NT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_WT_D26_DAC_vs_NT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/WT_D26_DAC_vs_NT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_WT_D26_BMN_vs_NT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/WT_D26_BMN_vs_NT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_WT_D26_DAC_BMN_vs_NT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/WT_D26_DAC_BMN_vs_NT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_KO_D5_DAC_vs_NT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/KO_D5_DAC_vs_NT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_KO_D5_BMN_vs_NT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/KO_D5_BMN_vs_NT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_KO_D5_DAC_BMN_vs_NT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/KO_D5_DAC_BMN_vs_NT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_KO_D16_DAC_vs_NT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/KO_D16_DAC_vs_NT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_KO_D16_BMN_vs_NT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/KO_D16_BMN_vs_NT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_KO_D16_DAC_BMN_vs_NT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/KO_D16_DAC_BMN_vs_NT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_KO_D26_DAC_vs_NT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/KO_D26_DAC_vs_NT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_KO_D26_BMN_vs_NT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/KO_D26_BMN_vs_NT_up.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_KO_D26_DAC_BMN_vs_NT_up,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/KO_D26_DAC_BMN_vs_NT_up.txt",row.names = T,col.names = T,quote = F)

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
write.table(filter_FDR_0.05_WT_D5_DAC_vs_NT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/WT_D5_DAC_vs_NT_down.txt",row.names = T,col.names = T,quote = F)
#write.table(filter_FDR_0.05_WT_D5_BMN_vs_NT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/WT_D5_BMN_vs_NT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_WT_D5_DAC_BMN_vs_NT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/WT_D5_DAC_BMN_vs_NT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_WT_D16_DAC_vs_NT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/WT_D16_DAC_vs_NT_down.txt",row.names = T,col.names = T,quote = F)
#write.table(filter_FDR_0.05_WT_D16_BMN_vs_NT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/WT_D16_BMN_vs_NT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_WT_D16_DAC_BMN_vs_NT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/WT_D16_DAC_BMN_vs_NT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_WT_D26_DAC_vs_NT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/WT_D26_DAC_vs_NT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_WT_D26_BMN_vs_NT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/WT_D26_BMN_vs_NT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_WT_D26_DAC_BMN_vs_NT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/WT_D26_DAC_BMN_vs_NT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_KO_D5_DAC_vs_NT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/KO_D5_DAC_vs_NT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_KO_D5_BMN_vs_NT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/KO_D5_BMN_vs_NT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_KO_D5_DAC_BMN_vs_NT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/KO_D5_DAC_BMN_vs_NT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_KO_D16_DAC_vs_NT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/KO_D16_DAC_vs_NT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_KO_D16_BMN_vs_NT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/KO_D16_BMN_vs_NT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_KO_D16_DAC_BMN_vs_NT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/KO_D16_DAC_BMN_vs_NT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_KO_D26_DAC_vs_NT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/KO_D26_DAC_vs_NT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_KO_D26_BMN_vs_NT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/KO_D26_BMN_vs_NT_down.txt",row.names = T,col.names = T,quote = F)
write.table(filter_FDR_0.05_KO_D26_DAC_BMN_vs_NT_down,"D:/Repeats_TE/transcriptid_redoTE/TEdata_output/KO_D26_DAC_BMN_vs_NT_down.txt",row.names = T,col.names = T,quote = F)
