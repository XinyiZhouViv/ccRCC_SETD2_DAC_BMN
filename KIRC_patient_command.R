#import patient and normal data, SETD2 status of patient data
data_mrna_seq_v2_rsem_normal_samples <- read.csv("data_mrna_seq_v2_rsem_normal_samples.csv",sep="\t", na.strings=c(""," "))
data_mrna_seq_v2_rsem_primary_tumor <- read.csv("data_mrna_seq_v2_rsem_primary_tumor.csv",sep="\t",na.strings=c(""," "))
PATIENT_SETD2_STATUS_SUMMARY <- read.csv("mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2).csv",sep="\t")

#devide patients into WT and MUtant
WT_PATIENTS <- PATIENT_SETD2_STATUS_SUMMARY[PATIENT_SETD2_STATUS_SUMMARY$Mutation=="WT",]
MU_PATIENTS <- PATIENT_SETD2_STATUS_SUMMARY[!PATIENT_SETD2_STATUS_SUMMARY$Mutation=="WT",]

#remove 2 NA SETD2 expressed patient data 
which(is.na(WT_PATIENTS), arr.ind=TRUE)
WT_PATIENTS <- na.omit(WT_PATIENTS)

#select WT patient ID based on their SETD2 expression level, pick out the least 25%
quantile(WT_PATIENTS$SETD2_mRNA_expression)
#      0%      25%      50%      75%     100% 
#184.429  963.480 1196.910 1420.780 2789.110
MID_HIGH_SETD2_PATIENT_ID <- WT_PATIENTS[WT_PATIENTS$SETD2_mRNA_expression>963.480,]
HIGH_SETD2_PATIENT_ID <- WT_PATIENTS[WT_PATIENTS$SETD2_mRNA_expression>1420.780,]
LOW_SETD2_PATIENT_ID <- rbind(WT_PATIENTS[WT_PATIENTS$SETD2_mRNA_expression<=963.480,],MU_PATIENTS)

#####create low SETD2 patient gene expression table, merge with normal tissue gene expression table
colnames(data_mrna_seq_v2_rsem_primary_tumor) %in% LOW_SETD2_PATIENT_ID$SAMPLE_ID
LOW_SETD2_PATIENT_gene_expression <- cbind(data_mrna_seq_v2_rsem_primary_tumor[,c(1,2)],data_mrna_seq_v2_rsem_primary_tumor[,colnames(data_mrna_seq_v2_rsem_primary_tumor) %in% LOW_SETD2_PATIENT_ID$SAMPLE_ID])

#merge patient and normal together, remove 13 NA and 37 duplications
x <- merge(LOW_SETD2_PATIENT_gene_expression,data_mrna_seq_v2_rsem_normal_samples[,2:74],by.x = "Entrez_Gene_Id", by.y = "Entrez_Gene_Id")
x <- na.omit(x)
which(duplicated(x$Entrez_Gene_Id),arr.ind=TRUE)
x <- x[!duplicated(x$Entrez_Gene_Id),]
row.names(x) <- x[,1]
count_x<-x[,3:321]

#distinguish groups, 1 = tumor, 2 = normal
group=factor(ifelse(as.numeric(substr(colnames(count_x),14,15))<10,'tumor','normal'))
table(group)

#Normalization and remove low expressed genes using EdegeR
DGE_count_x <- DGEList(counts=count_x,group=group, lib.size = colSums(count_x))
keep_cpm <- rowSums(cpm(DGE_count_x)>0.5) >= 20
DGE_count_x<-DGE_count_x[keep_cpm, ,keep.lib.sizes=FALSE] #gene number from 20490 to 16781
DGE_count_x <- calcNormFactors(DGE_count_x)
DGE_count_x$samples #check normalization factors
DGE_count_x$counts

#buid a list to show each sample's type (tumor/normal)
colData <- data.frame(Sample_Type=group)
row.names(colData)<-colnames(x[,3:321])

#build DESeq2() formula
y <- DESeqDataSetFromMatrix(countData=ceiling(DGE_count_x$counts),
                            colData = colData,
                            design = ~Sample_Type)
y <- DESeq(y)
resultsNames(y)
res <- results(y,contrast = c("Sample_Type","tumor","normal"))
resOrdered <- res[order(res$log2FoldChange,decreasing = T),]
resOrdered <- as.data.frame(resOrdered)

#Add Symbol, Ensemble ID
resOrdered <- transform(resOrdered,ENTREZID=row.names(resOrdered))
resOrdered_geneid <- bitr(resOrdered$ENTREZID,fromType="ENTREZID", toType=c("ENSEMBL","SYMBOL"), OrgDb="org.Hs.eg.db")

#merge show ENTREZID, ENSEMBL, SYMBOL
resOrdered <- merge(resOrdered,resOrdered_geneid,by.x="ENTREZID",by.y="ENTREZID")
resOrdered <- resOrdered[,c(1,8,9,2:7)]

#check and remove duplication
which(duplicated(resOrdered$ENTREZID), arr.ind=TRUE)
resOrdered <- filter(resOrdered,!duplicated(resOrdered$ENSEMBL))
resOrdered <- filter(resOrdered,!duplicated(resOrdered$SYMBOL))

#build up and down regulated gene list (SETD2 low expressed tumor vs. Normal tissue)
SETD2_LOW_T_vs_N_up_gene <- resOrdered[resOrdered$log2FoldChange>1 & resOrdered$padj<0.05,];row.names(SETD2_LOW_T_vs_N_up_gene) <- SETD2_LOW_T_vs_N_up_gene$ENSEMBL
SETD2_LOW_T_vs_N_down_gene <- resOrdered[resOrdered$log2FoldChange<(-1)  & resOrdered$padj<0.05,];row.names(SETD2_LOW_T_vs_N_down_gene) <- SETD2_LOW_T_vs_N_down_gene$ENSEMBL

write.table(SETD2_LOW_T_vs_N_up_gene,"SETD2_LOW_T_vs_N_up_gene.txt")
write.table(SETD2_LOW_T_vs_N_down_gene,"SETD2_LOW_T_vs_N_down_gene.txt")

#build up and down-regulated gene bar graph in SETD2 low tumor vs. Normal 
SETD2_LOW_T_vs_N_gene <- rbind(SETD2_LOW_T_vs_N_up_gene,SETD2_LOW_T_vs_N_down_gene)
ggplot(SETD2_LOW_T_vs_N_gene,aes(x=reorder(ENSEMBL,log2FoldChange), y=log2FoldChange,fill=log2FoldChange))+
  geom_bar(stat="identity",show.legend = T)+
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(hjust = 0.55,size=10),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6),
        plot.title = element_text(hjust = 0.5,size = 12))+
  labs(x="",y="logFC",title = "Up & down regulated genes in low SETD2 patients (FC>2)")+
  scale_fill_gradient2(limits=c(-6,6),midpoint = 0,low = "royalblue",high = "red")

#flip
SETD2_LOW_T_vs_N_gene <- rbind(SETD2_LOW_T_vs_N_up_gene,SETD2_LOW_T_vs_N_down_gene)
ggplot(SETD2_LOW_T_vs_N_gene,aes(x=reorder(ENSEMBL,log2FoldChange), y=log2FoldChange,fill=log2FoldChange))+
  geom_bar(stat="identity",show.legend = T)+
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(hjust = 0.55,size=10),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6),
        plot.title = element_text(hjust = 0.5,size = 12))+
  labs(x="Up & down regulated genes in low SETD2 patients (FC>2)",y="logFC",title = "")+
  scale_fill_gradient2(limits=c(-6,6),midpoint = 0,low = "royalblue",high = "red")+
  coord_flip()

#####create HIGH SETD2 patient gene expression table, merge with normal tissue gene expression table
colnames(data_mrna_seq_v2_rsem_primary_tumor) %in% HIGH_SETD2_PATIENT_ID$SAMPLE_ID
HIGH_SETD2_PATIENT_gene_expression <- cbind(data_mrna_seq_v2_rsem_primary_tumor[,c(1,2)],data_mrna_seq_v2_rsem_primary_tumor[,colnames(data_mrna_seq_v2_rsem_primary_tumor) %in% HIGH_SETD2_PATIENT_ID$SAMPLE_ID])

#merge patient and normal together, remove 13 NA and 37 duplications
x_HIGH <- merge(HIGH_SETD2_PATIENT_gene_expression,data_mrna_seq_v2_rsem_normal_samples[,2:74],by.x = "Entrez_Gene_Id", by.y = "Entrez_Gene_Id")
x_HIGH <- na.omit(x_HIGH)
which(duplicated(x_HIGH$Entrez_Gene_Id),arr.ind=TRUE)
x_HIGH <- x_HIGH[!duplicated(x_HIGH$Entrez_Gene_Id),]
row.names(x_HIGH) <- x_HIGH[,1]
count_x_HIGH<-x_HIGH[,3:162]

#distinguish groups, 1 = tumor, 2 = normal
group_HIGH=factor(ifelse(as.numeric(substr(colnames(count_x_HIGH),14,15))<10,'tumor','normal'))
table(group_HIGH)

#Normalization and remove low expressed genes using EdegeR
DGE_count_x_HIGH <- DGEList(counts=count_x_HIGH,group=group_HIGH, lib.size = colSums(count_x_HIGH))
keep_cpm_HIGH <- rowSums(cpm(DGE_count_x_HIGH)>0.5) >= 10
DGE_count_x_HIGH<-DGE_count_x_HIGH[keep_cpm_HIGH, ,keep.lib.sizes=FALSE] #gene number from 20490 to 16734
DGE_count_x_HIGH <- calcNormFactors(DGE_count_x_HIGH)
DGE_count_x_HIGH$samples #check normalization factors
DGE_count_x_HIGH$counts

#buid a list to show each sample's type (tumor/normal)
colData_HIGH <- data.frame(Sample_Type=group_HIGH)
row.names(colData_HIGH)<-colnames(x_HIGH[,3:162])

#build DESeq2() formula
y_HIGH <- DESeqDataSetFromMatrix(countData=ceiling(DGE_count_x_HIGH$counts),
                            colData = colData_HIGH,
                            design = ~Sample_Type)
y_HIGH <- DESeq(y_HIGH)
resultsNames(y_HIGH)
res_HIGH <- results(y_HIGH,contrast = c("Sample_Type","tumor","normal"))
resOrdered_HIGH <- res_HIGH[order(res_HIGH$log2FoldChange,decreasing = T),]
resOrdered_HIGH <- as.data.frame(resOrdered_HIGH)

#Add Symbol, Ensemble ID
resOrdered_HIGH <- transform(resOrdered_HIGH,ENTREZID=row.names(resOrdered_HIGH))
resOrdered_HIGH_geneid <- bitr(resOrdered_HIGH$ENTREZID,fromType="ENTREZID", toType=c("ENSEMBL","SYMBOL"), OrgDb="org.Hs.eg.db")

#merge show ENTREZID, ENSEMBL, SYMBOL
resOrdered_HIGH <- merge(resOrdered_HIGH,resOrdered_HIGH_geneid,by.x="ENTREZID",by.y="ENTREZID")
resOrdered_HIGH <- resOrdered_HIGH[,c(1,8,9,2:7)]

#check and remove duplication
which(duplicated(resOrdered_HIGH$ENTREZID), arr.ind=TRUE)
resOrdered_HIGH <- filter(resOrdered_HIGH,!duplicated(resOrdered_HIGH$ENSEMBL))
resOrdered_HIGH <- filter(resOrdered_HIGH,!duplicated(resOrdered_HIGH$SYMBOL))#16362 GENES

#build up and down regulated gene list (SETD2 low expressed tumor vs. Normal tissue)
SETD2_HIGH_T_vs_N_up_gene <- resOrdered_HIGH[resOrdered_HIGH$log2FoldChange>1 & resOrdered_HIGH$padj<0.05,];row.names(SETD2_HIGH_T_vs_N_up_gene) <- SETD2_HIGH_T_vs_N_up_gene$ENSEMBL#2492 genes
SETD2_HIGH_T_vs_N_down_gene <- resOrdered_HIGH[resOrdered_HIGH$log2FoldChange<(-1)  & resOrdered_HIGH$padj<0.05,];row.names(SETD2_HIGH_T_vs_N_down_gene) <- SETD2_HIGH_T_vs_N_down_gene$ENSEMBL #1613 genes

write.table(SETD2_HIGH_T_vs_N_up_gene,"SETD2_HIGH_T_vs_N_up_gene.txt")
write.table(SETD2_HIGH_T_vs_N_down_gene,"SETD2_HIGH_T_vs_N_down_gene.txt")

#build up and down-regulated gene bar graph in SETD2 high tumor vs. Normal 
SETD2_HIGH_T_vs_N_gene <- rbind(SETD2_HIGH_T_vs_N_up_gene,SETD2_HIGH_T_vs_N_down_gene)
ggplot(SETD2_HIGH_T_vs_N_gene,aes(x=reorder(ENSEMBL,log2FoldChange), y=log2FoldChange,fill=log2FoldChange))+
  geom_bar(stat="identity",show.legend = T)+
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(hjust = 0.55,size=10),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6),
        plot.title = element_text(hjust = 0.5,size = 12))+
  labs(x="",y="logFC",title = "Up & down regulated genes in high SETD2 patients (FC>2)")+
  scale_fill_gradient2(limits=c(-6,6),midpoint = 0,low = "royalblue",high = "red")

#flip
SETD2_LOW_T_vs_N_gene <- rbind(SETD2_LOW_T_vs_N_up_gene,SETD2_LOW_T_vs_N_down_gene)
ggplot(SETD2_HIGH_T_vs_N_gene,aes(x=reorder(ENSEMBL,log2FoldChange), y=log2FoldChange,fill=log2FoldChange))+
  geom_bar(stat="identity",show.legend = T)+
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(hjust = 0.55,size=10),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6),
        plot.title = element_text(hjust = 0.5,size = 12))+
  labs(x="Up & down regulated genes in high SETD2 patients (FC>2)",y="logFC",title = "")+
  scale_fill_gradient2(limits=c(-6,6),midpoint = 0,low = "royalblue",high = "red")+
  coord_flip()
