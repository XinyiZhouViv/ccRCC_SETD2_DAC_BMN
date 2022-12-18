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


###dotplot_KIRC_patient_command
#import DEG result table under treatment
filtered_edgeRt_WT_D5_DAC_vs_NT <- filtered_edgeRt_WT_D5_DAC_vs_NT[filtered_edgeRt_WT_D5_DAC_vs_NT$FDR < 0.05,]
filtered_edgeRt_KO_D5_DAC_vs_NT <- filtered_edgeRt_KO_D5_DAC_vs_NT[filtered_edgeRt_KO_D5_DAC_vs_NT$FDR < 0.05,]
filtered_edgeRt_WT_D5_BMN_vs_NT <- filtered_edgeRt_WT_D5_BMN_vs_NT[filtered_edgeRt_WT_D5_BMN_vs_NT$FDR < 0.05,]
filtered_edgeRt_KO_D5_BMN_vs_NT <- filtered_edgeRt_KO_D5_BMN_vs_NT[filtered_edgeRt_KO_D5_BMN_vs_NT$FDR < 0.05,]
filtered_edgeRt_WT_D5_DAC_BMN_vs_NT <- filtered_edgeRt_WT_D5_DAC_BMN_vs_NT[filtered_edgeRt_WT_D5_DAC_BMN_vs_NT$FDR < 0.05,]
filtered_edgeRt_KO_D5_DAC_BMN_vs_NT <- filtered_edgeRt_KO_D5_DAC_BMN_vs_NT[filtered_edgeRt_KO_D5_DAC_BMN_vs_NT$FDR < 0.05,]
filtered_edgeRt_WT_D16_DAC_vs_NT <- filtered_edgeRt_WT_D16_DAC_vs_NT[filtered_edgeRt_WT_D16_DAC_vs_NT$FDR < 0.05,]
filtered_edgeRt_KO_D16_DAC_vs_NT <- filtered_edgeRt_KO_D16_DAC_vs_NT[filtered_edgeRt_KO_D16_DAC_vs_NT$FDR < 0.05,]
filtered_edgeRt_WT_D16_BMN_vs_NT <- filtered_edgeRt_WT_D16_BMN_vs_NT[filtered_edgeRt_WT_D16_BMN_vs_NT$FDR < 0.05,]
filtered_edgeRt_KO_D16_BMN_vs_NT <- filtered_edgeRt_KO_D16_BMN_vs_NT[filtered_edgeRt_KO_D16_BMN_vs_NT$FDR < 0.05,]
filtered_edgeRt_WT_D16_DAC_BMN_vs_NT <- filtered_edgeRt_WT_D16_DAC_BMN_vs_NT[filtered_edgeRt_WT_D16_DAC_BMN_vs_NT$FDR < 0.05,]
filtered_edgeRt_KO_D16_DAC_BMN_vs_NT <- filtered_edgeRt_KO_D16_DAC_BMN_vs_NT[filtered_edgeRt_KO_D16_DAC_BMN_vs_NT$FDR < 0.05,]
filtered_edgeRt_WT_D26_DAC_vs_NT <- filtered_edgeRt_WT_D26_DAC_vs_NT[filtered_edgeRt_WT_D26_DAC_vs_NT$FDR < 0.05,]
filtered_edgeRt_KO_D26_DAC_vs_NT <- filtered_edgeRt_KO_D26_DAC_vs_NT[filtered_edgeRt_KO_D26_DAC_vs_NT$FDR < 0.05,]
filtered_edgeRt_WT_D26_BMN_vs_NT <- filtered_edgeRt_WT_D26_BMN_vs_NT[filtered_edgeRt_WT_D26_BMN_vs_NT$FDR < 0.05,]
filtered_edgeRt_KO_D26_BMN_vs_NT <- filtered_edgeRt_KO_D26_BMN_vs_NT[filtered_edgeRt_KO_D26_BMN_vs_NT$FDR < 0.05,]
filtered_edgeRt_WT_D26_DAC_BMN_vs_NT <- filtered_edgeRt_WT_D26_DAC_BMN_vs_NT[filtered_edgeRt_WT_D26_DAC_BMN_vs_NT$FDR < 0.05,]
filtered_edgeRt_KO_D26_DAC_BMN_vs_NT <- filtered_edgeRt_KO_D26_DAC_BMN_vs_NT[filtered_edgeRt_KO_D26_DAC_BMN_vs_NT$FDR < 0.05,]

#screen out each comparison's logFC of genes in (SETD2 high/low patients) UP regulated gene list
#Day5
WT_DAC_D5_FC_inuplist <- filtered_edgeRt_WT_D5_DAC_vs_NT;WT_DAC_D5_FC_inuplist <- transform(WT_DAC_D5_FC_inuplist,inpatients=WT_DAC_D5_FC_inuplist$X%in%row.names(SETD2_HIGH_T_vs_N_up_gene));WT_DAC_D5_FC_inuplist <- WT_DAC_D5_FC_inuplist[WT_DAC_D5_FC_inuplist$inpatients==T,];WT_DAC_D5_FC_inuplist <- WT_DAC_D5_FC_inuplist[,c(1,3)];WT_DAC_D5_FC_inuplist <- transform(WT_DAC_D5_FC_inuplist,Treatment="DAC");WT_DAC_D5_FC_inuplist <- transform(WT_DAC_D5_FC_inuplist,Celltype="WT");WT_DAC_D5_FC_inuplist <- WT_DAC_D5_FC_inuplist[order(WT_DAC_D5_FC_inuplist$logFC,decreasing=T),]

WT_BMN_D5_FC_inuplist <- filtered_edgeRt_WT_D5_BMN_vs_NT;WT_BMN_D5_FC_inuplist <- transform(WT_BMN_D5_FC_inuplist,inpatients=WT_BMN_D5_FC_inuplist$X%in%row.names(SETD2_HIGH_T_vs_N_up_gene));WT_BMN_D5_FC_inuplist <- WT_BMN_D5_FC_inuplist[WT_BMN_D5_FC_inuplist$inpatients==T,];WT_BMN_D5_FC_inuplist <- WT_BMN_D5_FC_inuplist[,c(1,3)];WT_BMN_D5_FC_inuplist <- transform(WT_BMN_D5_FC_inuplist,Treatment="BMN");WT_BMN_D5_FC_inuplist <- transform(WT_BMN_D5_FC_inuplist,Celltype="WT")

WT_DAC_BMN_D5_FC_inuplist <- filtered_edgeRt_WT_D5_DAC_BMN_vs_NT;WT_DAC_BMN_D5_FC_inuplist <- transform(WT_DAC_BMN_D5_FC_inuplist,inpatients=WT_DAC_BMN_D5_FC_inuplist$X%in%row.names(SETD2_HIGH_T_vs_N_up_gene));WT_DAC_BMN_D5_FC_inuplist <- WT_DAC_BMN_D5_FC_inuplist[WT_DAC_BMN_D5_FC_inuplist$inpatients==T,];WT_DAC_BMN_D5_FC_inuplist <- WT_DAC_BMN_D5_FC_inuplist[,c(1,3)];WT_DAC_BMN_D5_FC_inuplist <- transform(WT_DAC_BMN_D5_FC_inuplist,Treatment="DAC_BMN");WT_DAC_BMN_D5_FC_inuplist <- transform(WT_DAC_BMN_D5_FC_inuplist,Celltype="WT")

KO_DAC_D5_FC_inuplist <- filtered_edgeRt_KO_D5_DAC_vs_NT;KO_DAC_D5_FC_inuplist <- transform(KO_DAC_D5_FC_inuplist,inpatients=KO_DAC_D5_FC_inuplist$X%in%row.names(SETD2_LOW_T_vs_N_up_gene));KO_DAC_D5_FC_inuplist <- KO_DAC_D5_FC_inuplist[KO_DAC_D5_FC_inuplist$inpatients==T,];KO_DAC_D5_FC_inuplist <- KO_DAC_D5_FC_inuplist[,c(1,3)];KO_DAC_D5_FC_inuplist <- transform(KO_DAC_D5_FC_inuplist,Treatment="DAC");KO_DAC_D5_FC_inuplist <- transform(KO_DAC_D5_FC_inuplist,Celltype="KO");KO_DAC_D5_FC_inuplist <- KO_DAC_D5_FC_inuplist[order(KO_DAC_D5_FC_inuplist$logFC,decreasing=T),]

KO_BMN_D5_FC_inuplist <- filtered_edgeRt_KO_D5_BMN_vs_NT;KO_BMN_D5_FC_inuplist <- transform(KO_BMN_D5_FC_inuplist,inpatients=KO_BMN_D5_FC_inuplist$X%in%row.names(SETD2_LOW_T_vs_N_up_gene));KO_BMN_D5_FC_inuplist <- KO_BMN_D5_FC_inuplist[KO_BMN_D5_FC_inuplist$inpatients==T,];KO_BMN_D5_FC_inuplist <- KO_BMN_D5_FC_inuplist[,c(1,3)];KO_BMN_D5_FC_inuplist <- transform(KO_BMN_D5_FC_inuplist,Treatment="BMN");KO_BMN_D5_FC_inuplist <- transform(KO_BMN_D5_FC_inuplist,Celltype="KO")

KO_DAC_BMN_D5_FC_inuplist <- filtered_edgeRt_KO_D5_DAC_BMN_vs_NT;KO_DAC_BMN_D5_FC_inuplist <- transform(KO_DAC_BMN_D5_FC_inuplist,inpatients=KO_DAC_BMN_D5_FC_inuplist$X%in%row.names(SETD2_LOW_T_vs_N_up_gene));KO_DAC_BMN_D5_FC_inuplist <- KO_DAC_BMN_D5_FC_inuplist[KO_DAC_BMN_D5_FC_inuplist$inpatients==T,];KO_DAC_BMN_D5_FC_inuplist <- KO_DAC_BMN_D5_FC_inuplist[,c(1,3)];KO_DAC_BMN_D5_FC_inuplist <- transform(KO_DAC_BMN_D5_FC_inuplist,Treatment="DAC_BMN");KO_DAC_BMN_D5_FC_inuplist <- transform(KO_DAC_BMN_D5_FC_inuplist,Celltype="KO")

# WT
Sequence_WT_D5_inuplist <- transform(WT_DAC_D5_FC_inuplist,Genesequence=c(1:644));Sequence_WT_D5_inuplist <- Sequence_WT_D5_inuplist[,c(1,5)]

Geneset_WT_D5_inuplist <- rbind(WT_DAC_D5_FC_inuplist,WT_BMN_D5_FC_inuplist,WT_DAC_BMN_D5_FC_inuplist)

merge_WT_up_D5 <- merge(Geneset_WT_D5_inuplist,Sequence_WT_D5_inuplist,by.x = "X",by.y = "X")

ggplot(merge_WT_up_D5,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes up-regulated in HIGH-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression change under treatments in WT 786-O cells (Day5)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank())+
  scale_color_manual(values = c("red","blue"))

# KO
Sequence_KO_D5_inuplist <- transform(KO_DAC_D5_FC_inuplist,Genesequence=c(1:1070));Sequence_KO_D5_inuplist <- Sequence_KO_D5_inuplist[,c(1,5)]

Geneset_KO_D5_inuplist <- rbind(KO_DAC_D5_FC_inuplist,KO_BMN_D5_FC_inuplist,KO_DAC_BMN_D5_FC_inuplist)

merge_KO_up_D5 <- merge(Geneset_KO_D5_inuplist,Sequence_KO_D5_inuplist,by.x = "X",by.y = "X")

ggplot(merge_KO_up_D5,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes up-regulated in LOW-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression change under treatments in KO 786-O cells (Day5)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank())+
scale_color_manual(values = c("green","red","blue"))

#Day16
WT_DAC_D16_FC_inuplist <- filtered_edgeRt_WT_D16_DAC_vs_NT;WT_DAC_D16_FC_inuplist <- transform(WT_DAC_D16_FC_inuplist,inpatients=WT_DAC_D16_FC_inuplist$X%in%row.names(SETD2_HIGH_T_vs_N_up_gene));WT_DAC_D16_FC_inuplist <- WT_DAC_D16_FC_inuplist[WT_DAC_D16_FC_inuplist$inpatients==T,];WT_DAC_D16_FC_inuplist <- WT_DAC_D16_FC_inuplist[,c(1,3)];WT_DAC_D16_FC_inuplist <- transform(WT_DAC_D16_FC_inuplist,Treatment="DAC");WT_DAC_D16_FC_inuplist <- transform(WT_DAC_D16_FC_inuplist,Celltype="WT");WT_DAC_D16_FC_inuplist <- WT_DAC_D16_FC_inuplist[order(WT_DAC_D16_FC_inuplist$logFC,decreasing=T),]

WT_BMN_D16_FC_inuplist <- filtered_edgeRt_WT_D16_BMN_vs_NT;WT_BMN_D16_FC_inuplist <- transform(WT_BMN_D16_FC_inuplist,inpatients=WT_BMN_D16_FC_inuplist$X%in%row.names(SETD2_HIGH_T_vs_N_up_gene));WT_BMN_D16_FC_inuplist <- WT_BMN_D16_FC_inuplist[WT_BMN_D16_FC_inuplist$inpatients==T,];WT_BMN_D16_FC_inuplist <- WT_BMN_D16_FC_inuplist[,c(1,3)];WT_BMN_D16_FC_inuplist <- transform(WT_BMN_D16_FC_inuplist,Treatment="BMN");WT_BMN_D16_FC_inuplist <- transform(WT_BMN_D16_FC_inuplist,Celltype="WT")

WT_DAC_BMN_D16_FC_inuplist <- filtered_edgeRt_WT_D16_DAC_BMN_vs_NT;WT_DAC_BMN_D16_FC_inuplist <- transform(WT_DAC_BMN_D16_FC_inuplist,inpatients=WT_DAC_BMN_D16_FC_inuplist$X%in%row.names(SETD2_HIGH_T_vs_N_up_gene));WT_DAC_BMN_D16_FC_inuplist <- WT_DAC_BMN_D16_FC_inuplist[WT_DAC_BMN_D16_FC_inuplist$inpatients==T,];WT_DAC_BMN_D16_FC_inuplist <- WT_DAC_BMN_D16_FC_inuplist[,c(1,3)];WT_DAC_BMN_D16_FC_inuplist <- transform(WT_DAC_BMN_D16_FC_inuplist,Treatment="DAC_BMN");WT_DAC_BMN_D16_FC_inuplist <- transform(WT_DAC_BMN_D16_FC_inuplist,Celltype="WT")

KO_DAC_D16_FC_inuplist <- filtered_edgeRt_KO_D16_DAC_vs_NT;KO_DAC_D16_FC_inuplist <- transform(KO_DAC_D16_FC_inuplist,inpatients=KO_DAC_D16_FC_inuplist$X%in%row.names(SETD2_LOW_T_vs_N_up_gene));KO_DAC_D16_FC_inuplist <- KO_DAC_D16_FC_inuplist[KO_DAC_D16_FC_inuplist$inpatients==T,];KO_DAC_D16_FC_inuplist <- KO_DAC_D16_FC_inuplist[,c(1,3)];KO_DAC_D16_FC_inuplist <- transform(KO_DAC_D16_FC_inuplist,Treatment="DAC");KO_DAC_D16_FC_inuplist <- transform(KO_DAC_D16_FC_inuplist,Celltype="KO");KO_DAC_D16_FC_inuplist <- KO_DAC_D16_FC_inuplist[order(KO_DAC_D16_FC_inuplist$logFC,decreasing=T),]

KO_BMN_D16_FC_inuplist <- filtered_edgeRt_KO_D16_BMN_vs_NT;KO_BMN_D16_FC_inuplist <- transform(KO_BMN_D16_FC_inuplist,inpatients=KO_BMN_D16_FC_inuplist$X%in%row.names(SETD2_LOW_T_vs_N_up_gene));KO_BMN_D16_FC_inuplist <- KO_BMN_D16_FC_inuplist[KO_BMN_D16_FC_inuplist$inpatients==T,];KO_BMN_D16_FC_inuplist <- KO_BMN_D16_FC_inuplist[,c(1,3)];KO_BMN_D16_FC_inuplist <- transform(KO_BMN_D16_FC_inuplist,Treatment="BMN");KO_BMN_D16_FC_inuplist <- transform(KO_BMN_D16_FC_inuplist,Celltype="KO")

KO_DAC_BMN_D16_FC_inuplist <- filtered_edgeRt_KO_D16_DAC_BMN_vs_NT;KO_DAC_BMN_D16_FC_inuplist <- transform(KO_DAC_BMN_D16_FC_inuplist,inpatients=KO_DAC_BMN_D16_FC_inuplist$X%in%row.names(SETD2_LOW_T_vs_N_up_gene));KO_DAC_BMN_D16_FC_inuplist <- KO_DAC_BMN_D16_FC_inuplist[KO_DAC_BMN_D16_FC_inuplist$inpatients==T,];KO_DAC_BMN_D16_FC_inuplist <- KO_DAC_BMN_D16_FC_inuplist[,c(1,3)];KO_DAC_BMN_D16_FC_inuplist <- transform(KO_DAC_BMN_D16_FC_inuplist,Treatment="DAC_BMN");KO_DAC_BMN_D16_FC_inuplist <- transform(KO_DAC_BMN_D16_FC_inuplist,Celltype="KO")

# WT
Sequence_WT_D16_inuplist <- transform(WT_DAC_D16_FC_inuplist,Genesequence=c(1:881));Sequence_WT_D16_inuplist <- Sequence_WT_D16_inuplist[,c(1,5)]

Geneset_WT_D16_inuplist <- rbind(WT_DAC_D16_FC_inuplist,WT_BMN_D16_FC_inuplist,WT_DAC_BMN_D16_FC_inuplist)

merge_WT_up_D16 <- merge(Geneset_WT_D16_inuplist,Sequence_WT_D16_inuplist,by.x = "X",by.y = "X")

ggplot(merge_WT_up_D16,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes up-regulated in HIGH-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression change under treatments in WT 786-O cells (Day16)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank())+
  scale_color_manual(values = c("green","red","blue"))

# KO
Sequence_KO_D16_inuplist <- transform(KO_DAC_D16_FC_inuplist,Genesequence=c(1:1206));Sequence_KO_D16_inuplist <- Sequence_KO_D16_inuplist[,c(1,5)]

Geneset_KO_D16_inuplist <- rbind(KO_DAC_D16_FC_inuplist,KO_BMN_D16_FC_inuplist,KO_DAC_BMN_D16_FC_inuplist)

merge_KO_up_D16 <- merge(Geneset_KO_D16_inuplist,Sequence_KO_D16_inuplist,by.x = "X",by.y = "X")

ggplot(merge_KO_up_D16,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes up-regulated in LOW-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression change under treatments in KO 786-O cells (Day16)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank())+
  scale_color_manual(values = c("green","red","blue"))


#Day26
WT_DAC_D26_FC_inuplist <- filtered_edgeRt_WT_D26_DAC_vs_NT;WT_DAC_D26_FC_inuplist <- transform(WT_DAC_D26_FC_inuplist,inpatients=WT_DAC_D26_FC_inuplist$X%in%row.names(SETD2_HIGH_T_vs_N_up_gene));WT_DAC_D26_FC_inuplist <- WT_DAC_D26_FC_inuplist[WT_DAC_D26_FC_inuplist$inpatients==T,];WT_DAC_D26_FC_inuplist <- WT_DAC_D26_FC_inuplist[,c(1,3)];WT_DAC_D26_FC_inuplist <- transform(WT_DAC_D26_FC_inuplist,Treatment="DAC");WT_DAC_D26_FC_inuplist <- transform(WT_DAC_D26_FC_inuplist,Celltype="WT");WT_DAC_D26_FC_inuplist <- WT_DAC_D26_FC_inuplist[order(WT_DAC_D26_FC_inuplist$logFC,decreasing=T),]

WT_BMN_D26_FC_inuplist <- filtered_edgeRt_WT_D26_BMN_vs_NT;WT_BMN_D26_FC_inuplist <- transform(WT_BMN_D26_FC_inuplist,inpatients=WT_BMN_D26_FC_inuplist$X%in%row.names(SETD2_HIGH_T_vs_N_up_gene));WT_BMN_D26_FC_inuplist <- WT_BMN_D26_FC_inuplist[WT_BMN_D26_FC_inuplist$inpatients==T,];WT_BMN_D26_FC_inuplist <- WT_BMN_D26_FC_inuplist[,c(1,3)];WT_BMN_D26_FC_inuplist <- transform(WT_BMN_D26_FC_inuplist,Treatment="BMN");WT_BMN_D26_FC_inuplist <- transform(WT_BMN_D26_FC_inuplist,Celltype="WT")

WT_DAC_BMN_D26_FC_inuplist <- filtered_edgeRt_WT_D26_DAC_BMN_vs_NT;WT_DAC_BMN_D26_FC_inuplist <- transform(WT_DAC_BMN_D26_FC_inuplist,inpatients=WT_DAC_BMN_D26_FC_inuplist$X%in%row.names(SETD2_HIGH_T_vs_N_up_gene));WT_DAC_BMN_D26_FC_inuplist <- WT_DAC_BMN_D26_FC_inuplist[WT_DAC_BMN_D26_FC_inuplist$inpatients==T,];WT_DAC_BMN_D26_FC_inuplist <- WT_DAC_BMN_D26_FC_inuplist[,c(1,3)];WT_DAC_BMN_D26_FC_inuplist <- transform(WT_DAC_BMN_D26_FC_inuplist,Treatment="DAC_BMN");WT_DAC_BMN_D26_FC_inuplist <- transform(WT_DAC_BMN_D26_FC_inuplist,Celltype="WT")

KO_DAC_D26_FC_inuplist <- filtered_edgeRt_KO_D26_DAC_vs_NT;KO_DAC_D26_FC_inuplist <- transform(KO_DAC_D26_FC_inuplist,inpatients=KO_DAC_D26_FC_inuplist$X%in%row.names(SETD2_LOW_T_vs_N_up_gene));KO_DAC_D26_FC_inuplist <- KO_DAC_D26_FC_inuplist[KO_DAC_D26_FC_inuplist$inpatients==T,];KO_DAC_D26_FC_inuplist <- KO_DAC_D26_FC_inuplist[,c(1,3)];KO_DAC_D26_FC_inuplist <- transform(KO_DAC_D26_FC_inuplist,Treatment="DAC");KO_DAC_D26_FC_inuplist <- transform(KO_DAC_D26_FC_inuplist,Celltype="KO");KO_DAC_D26_FC_inuplist <- KO_DAC_D26_FC_inuplist[order(KO_DAC_D26_FC_inuplist$logFC,decreasing=T),]

KO_BMN_D26_FC_inuplist <- filtered_edgeRt_KO_D26_BMN_vs_NT;KO_BMN_D26_FC_inuplist <- transform(KO_BMN_D26_FC_inuplist,inpatients=KO_BMN_D26_FC_inuplist$X%in%row.names(SETD2_LOW_T_vs_N_up_gene));KO_BMN_D26_FC_inuplist <- KO_BMN_D26_FC_inuplist[KO_BMN_D26_FC_inuplist$inpatients==T,];KO_BMN_D26_FC_inuplist <- KO_BMN_D26_FC_inuplist[,c(1,3)];KO_BMN_D26_FC_inuplist <- transform(KO_BMN_D26_FC_inuplist,Treatment="BMN");KO_BMN_D26_FC_inuplist <- transform(KO_BMN_D26_FC_inuplist,Celltype="KO")

KO_DAC_BMN_D26_FC_inuplist <- filtered_edgeRt_KO_D26_DAC_BMN_vs_NT;KO_DAC_BMN_D26_FC_inuplist <- transform(KO_DAC_BMN_D26_FC_inuplist,inpatients=KO_DAC_BMN_D26_FC_inuplist$X%in%row.names(SETD2_LOW_T_vs_N_up_gene));KO_DAC_BMN_D26_FC_inuplist <- KO_DAC_BMN_D26_FC_inuplist[KO_DAC_BMN_D26_FC_inuplist$inpatients==T,];KO_DAC_BMN_D26_FC_inuplist <- KO_DAC_BMN_D26_FC_inuplist[,c(1,3)];KO_DAC_BMN_D26_FC_inuplist <- transform(KO_DAC_BMN_D26_FC_inuplist,Treatment="DAC_BMN");KO_DAC_BMN_D26_FC_inuplist <- transform(KO_DAC_BMN_D26_FC_inuplist,Celltype="KO")

# WT
Sequence_WT_D26_inuplist <- transform(WT_DAC_D26_FC_inuplist,Genesequence=c(1:643));Sequence_WT_D26_inuplist <- Sequence_WT_D26_inuplist[,c(1,5)]

Geneset_WT_D26_inuplist <- rbind(WT_DAC_D26_FC_inuplist,WT_BMN_D26_FC_inuplist,WT_DAC_BMN_D26_FC_inuplist)

merge_WT_up_D26 <- merge(Geneset_WT_D26_inuplist,Sequence_WT_D26_inuplist,by.x = "X",by.y = "X")

ggplot(merge_WT_up_D26,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes up-regulated in HIGH-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression change under treatments in WT 786-O cells (Day26)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank())+
  scale_color_manual(values = c("green","red","blue"))

# KO
Sequence_KO_D26_inuplist <- transform(KO_DAC_D26_FC_inuplist,Genesequence=c(1:1272));Sequence_KO_D26_inuplist <- Sequence_KO_D26_inuplist[,c(1,5)]

Geneset_KO_D26_inuplist <- rbind(KO_DAC_D26_FC_inuplist,KO_BMN_D26_FC_inuplist,KO_DAC_BMN_D26_FC_inuplist)

merge_KO_up_D26 <- merge(Geneset_KO_D26_inuplist,Sequence_KO_D26_inuplist,by.x = "X",by.y = "X")

ggplot(merge_KO_up_D26,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes up-regulated in LOW-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression change under treatments in KO 786-O cells (Day26)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank())+
  scale_color_manual(values = c("green","red","blue"))


#screen out each comparison's logFC of genes in (SETD2 high/low patients) DOWN regulated gene list
#Day5
WT_DAC_D5_FC_indownlist <- filtered_edgeRt_WT_D5_DAC_vs_NT;WT_DAC_D5_FC_indownlist <- transform(WT_DAC_D5_FC_indownlist,inpatients=WT_DAC_D5_FC_indownlist$X%in%row.names(SETD2_HIGH_T_vs_N_down_gene));WT_DAC_D5_FC_indownlist <- WT_DAC_D5_FC_indownlist[WT_DAC_D5_FC_indownlist$inpatients==T,];WT_DAC_D5_FC_indownlist <- WT_DAC_D5_FC_indownlist[,c(1,3)];WT_DAC_D5_FC_indownlist <- transform(WT_DAC_D5_FC_indownlist,Treatment="DAC");WT_DAC_D5_FC_indownlist <- transform(WT_DAC_D5_FC_indownlist,Celltype="WT");WT_DAC_D5_FC_indownlist <- WT_DAC_D5_FC_indownlist[order(WT_DAC_D5_FC_indownlist$logFC,decreasing=T),]

WT_BMN_D5_FC_indownlist <- filtered_edgeRt_WT_D5_BMN_vs_NT;WT_BMN_D5_FC_indownlist <- transform(WT_BMN_D5_FC_indownlist,inpatients=WT_BMN_D5_FC_indownlist$X%in%row.names(SETD2_HIGH_T_vs_N_down_gene));WT_BMN_D5_FC_indownlist <- WT_BMN_D5_FC_indownlist[WT_BMN_D5_FC_indownlist$inpatients==T,];WT_BMN_D5_FC_indownlist <- WT_BMN_D5_FC_indownlist[,c(1,3)];WT_BMN_D5_FC_indownlist <- transform(WT_BMN_D5_FC_indownlist,Treatment="BMN");WT_BMN_D5_FC_indownlist <- transform(WT_BMN_D5_FC_indownlist,Celltype="WT")

WT_DAC_BMN_D5_FC_indownlist <- filtered_edgeRt_WT_D5_DAC_BMN_vs_NT;WT_DAC_BMN_D5_FC_indownlist <- transform(WT_DAC_BMN_D5_FC_indownlist,inpatients=WT_DAC_BMN_D5_FC_indownlist$X%in%row.names(SETD2_HIGH_T_vs_N_down_gene));WT_DAC_BMN_D5_FC_indownlist <- WT_DAC_BMN_D5_FC_indownlist[WT_DAC_BMN_D5_FC_indownlist$inpatients==T,];WT_DAC_BMN_D5_FC_indownlist <- WT_DAC_BMN_D5_FC_indownlist[,c(1,3)];WT_DAC_BMN_D5_FC_indownlist <- transform(WT_DAC_BMN_D5_FC_indownlist,Treatment="DAC_BMN");WT_DAC_BMN_D5_FC_indownlist <- transform(WT_DAC_BMN_D5_FC_indownlist,Celltype="WT")

KO_DAC_D5_FC_indownlist <- filtered_edgeRt_KO_D5_DAC_vs_NT;KO_DAC_D5_FC_indownlist <- transform(KO_DAC_D5_FC_indownlist,inpatients=KO_DAC_D5_FC_indownlist$X%in%row.names(SETD2_LOW_T_vs_N_down_gene));KO_DAC_D5_FC_indownlist <- KO_DAC_D5_FC_indownlist[KO_DAC_D5_FC_indownlist$inpatients==T,];KO_DAC_D5_FC_indownlist <- KO_DAC_D5_FC_indownlist[,c(1,3)];KO_DAC_D5_FC_indownlist <- transform(KO_DAC_D5_FC_indownlist,Treatment="DAC");KO_DAC_D5_FC_indownlist <- transform(KO_DAC_D5_FC_indownlist,Celltype="KO");KO_DAC_D5_FC_indownlist <- KO_DAC_D5_FC_indownlist[order(KO_DAC_D5_FC_indownlist$logFC,decreasing=T),]

KO_BMN_D5_FC_indownlist <- filtered_edgeRt_KO_D5_BMN_vs_NT;KO_BMN_D5_FC_indownlist <- transform(KO_BMN_D5_FC_indownlist,inpatients=KO_BMN_D5_FC_indownlist$X%in%row.names(SETD2_LOW_T_vs_N_down_gene));KO_BMN_D5_FC_indownlist <- KO_BMN_D5_FC_indownlist[KO_BMN_D5_FC_indownlist$inpatients==T,];KO_BMN_D5_FC_indownlist <- KO_BMN_D5_FC_indownlist[,c(1,3)];KO_BMN_D5_FC_indownlist <- transform(KO_BMN_D5_FC_indownlist,Treatment="BMN");KO_BMN_D5_FC_indownlist <- transform(KO_BMN_D5_FC_indownlist,Celltype="KO")

KO_DAC_BMN_D5_FC_indownlist <- filtered_edgeRt_KO_D5_DAC_BMN_vs_NT;KO_DAC_BMN_D5_FC_indownlist <- transform(KO_DAC_BMN_D5_FC_indownlist,inpatients=KO_DAC_BMN_D5_FC_indownlist$X%in%row.names(SETD2_LOW_T_vs_N_down_gene));KO_DAC_BMN_D5_FC_indownlist <- KO_DAC_BMN_D5_FC_indownlist[KO_DAC_BMN_D5_FC_indownlist$inpatients==T,];KO_DAC_BMN_D5_FC_indownlist <- KO_DAC_BMN_D5_FC_indownlist[,c(1,3)];KO_DAC_BMN_D5_FC_indownlist <- transform(KO_DAC_BMN_D5_FC_indownlist,Treatment="DAC_BMN");KO_DAC_BMN_D5_FC_indownlist <- transform(KO_DAC_BMN_D5_FC_indownlist,Celltype="KO")

# WT
Sequence_WT_D5_indownlist <- transform(WT_DAC_D5_FC_indownlist,Genesequence=c(1:352));Sequence_WT_D5_indownlist <- Sequence_WT_D5_indownlist[,c(1,5)]

Geneset_WT_D5_indownlist <- rbind(WT_DAC_D5_FC_indownlist,WT_BMN_D5_FC_indownlist,WT_DAC_BMN_D5_FC_indownlist)

merge_WT_down_D5 <- merge(Geneset_WT_D5_indownlist,Sequence_WT_D5_indownlist,by.x = "X",by.y = "X")

ggplot(merge_WT_down_D5,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes down-regulated in HIGH-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression change under treatments in WT 786-O cells (Day5)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank())+
  scale_color_manual(values = c("red","blue"))

# KO
Sequence_KO_D5_indownlist <- transform(KO_DAC_D5_FC_indownlist,Genesequence=c(1:597));Sequence_KO_D5_indownlist <- Sequence_KO_D5_indownlist[,c(1,5)]

Geneset_KO_D5_indownlist <- rbind(KO_DAC_D5_FC_indownlist,KO_BMN_D5_FC_indownlist,KO_DAC_BMN_D5_FC_indownlist)

merge_KO_down_D5 <- merge(Geneset_KO_D5_indownlist,Sequence_KO_D5_indownlist,by.x = "X",by.y = "X")

ggplot(merge_KO_down_D5,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes down-regulated in LOW-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression change under treatments in KO 786-O cells (Day5)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank())+
  scale_color_manual(values = c("green","red","blue"))

#Day16
WT_DAC_D16_FC_indownlist <- filtered_edgeRt_WT_D16_DAC_vs_NT;WT_DAC_D16_FC_indownlist <- transform(WT_DAC_D16_FC_indownlist,inpatients=WT_DAC_D16_FC_indownlist$X%in%row.names(SETD2_HIGH_T_vs_N_down_gene));WT_DAC_D16_FC_indownlist <- WT_DAC_D16_FC_indownlist[WT_DAC_D16_FC_indownlist$inpatients==T,];WT_DAC_D16_FC_indownlist <- WT_DAC_D16_FC_indownlist[,c(1,3)];WT_DAC_D16_FC_indownlist <- transform(WT_DAC_D16_FC_indownlist,Treatment="DAC");WT_DAC_D16_FC_indownlist <- transform(WT_DAC_D16_FC_indownlist,Celltype="WT");WT_DAC_D16_FC_indownlist <- WT_DAC_D16_FC_indownlist[order(WT_DAC_D16_FC_indownlist$logFC,decreasing=T),]

WT_BMN_D16_FC_indownlist <- filtered_edgeRt_WT_D16_BMN_vs_NT;WT_BMN_D16_FC_indownlist <- transform(WT_BMN_D16_FC_indownlist,inpatients=WT_BMN_D16_FC_indownlist$X%in%row.names(SETD2_HIGH_T_vs_N_down_gene));WT_BMN_D16_FC_indownlist <- WT_BMN_D16_FC_indownlist[WT_BMN_D16_FC_indownlist$inpatients==T,];WT_BMN_D16_FC_indownlist <- WT_BMN_D16_FC_indownlist[,c(1,3)];WT_BMN_D16_FC_indownlist <- transform(WT_BMN_D16_FC_indownlist,Treatment="BMN");WT_BMN_D16_FC_indownlist <- transform(WT_BMN_D16_FC_indownlist,Celltype="WT")

WT_DAC_BMN_D16_FC_indownlist <- filtered_edgeRt_WT_D16_DAC_BMN_vs_NT;WT_DAC_BMN_D16_FC_indownlist <- transform(WT_DAC_BMN_D16_FC_indownlist,inpatients=WT_DAC_BMN_D16_FC_indownlist$X%in%row.names(SETD2_HIGH_T_vs_N_down_gene));WT_DAC_BMN_D16_FC_indownlist <- WT_DAC_BMN_D16_FC_indownlist[WT_DAC_BMN_D16_FC_indownlist$inpatients==T,];WT_DAC_BMN_D16_FC_indownlist <- WT_DAC_BMN_D16_FC_indownlist[,c(1,3)];WT_DAC_BMN_D16_FC_indownlist <- transform(WT_DAC_BMN_D16_FC_indownlist,Treatment="DAC_BMN");WT_DAC_BMN_D16_FC_indownlist <- transform(WT_DAC_BMN_D16_FC_indownlist,Celltype="WT")

KO_DAC_D16_FC_indownlist <- filtered_edgeRt_KO_D16_DAC_vs_NT;KO_DAC_D16_FC_indownlist <- transform(KO_DAC_D16_FC_indownlist,inpatients=KO_DAC_D16_FC_indownlist$X%in%row.names(SETD2_LOW_T_vs_N_down_gene));KO_DAC_D16_FC_indownlist <- KO_DAC_D16_FC_indownlist[KO_DAC_D16_FC_indownlist$inpatients==T,];KO_DAC_D16_FC_indownlist <- KO_DAC_D16_FC_indownlist[,c(1,3)];KO_DAC_D16_FC_indownlist <- transform(KO_DAC_D16_FC_indownlist,Treatment="DAC");KO_DAC_D16_FC_indownlist <- transform(KO_DAC_D16_FC_indownlist,Celltype="KO");KO_DAC_D16_FC_indownlist <- KO_DAC_D16_FC_indownlist[order(KO_DAC_D16_FC_indownlist$logFC,decreasing=T),]

KO_BMN_D16_FC_indownlist <- filtered_edgeRt_KO_D16_BMN_vs_NT;KO_BMN_D16_FC_indownlist <- transform(KO_BMN_D16_FC_indownlist,inpatients=KO_BMN_D16_FC_indownlist$X%in%row.names(SETD2_LOW_T_vs_N_down_gene));KO_BMN_D16_FC_indownlist <- KO_BMN_D16_FC_indownlist[KO_BMN_D16_FC_indownlist$inpatients==T,];KO_BMN_D16_FC_indownlist <- KO_BMN_D16_FC_indownlist[,c(1,3)];KO_BMN_D16_FC_indownlist <- transform(KO_BMN_D16_FC_indownlist,Treatment="BMN");KO_BMN_D16_FC_indownlist <- transform(KO_BMN_D16_FC_indownlist,Celltype="KO")

KO_DAC_BMN_D16_FC_indownlist <- filtered_edgeRt_KO_D16_DAC_BMN_vs_NT;KO_DAC_BMN_D16_FC_indownlist <- transform(KO_DAC_BMN_D16_FC_indownlist,inpatients=KO_DAC_BMN_D16_FC_indownlist$X%in%row.names(SETD2_LOW_T_vs_N_down_gene));KO_DAC_BMN_D16_FC_indownlist <- KO_DAC_BMN_D16_FC_indownlist[KO_DAC_BMN_D16_FC_indownlist$inpatients==T,];KO_DAC_BMN_D16_FC_indownlist <- KO_DAC_BMN_D16_FC_indownlist[,c(1,3)];KO_DAC_BMN_D16_FC_indownlist <- transform(KO_DAC_BMN_D16_FC_indownlist,Treatment="DAC_BMN");KO_DAC_BMN_D16_FC_indownlist <- transform(KO_DAC_BMN_D16_FC_indownlist,Celltype="KO")

# WT
Sequence_WT_D16_indownlist <- transform(WT_DAC_D16_FC_indownlist,Genesequence=c(1:423));Sequence_WT_D16_indownlist <- Sequence_WT_D16_indownlist[,c(1,5)]

Geneset_WT_D16_indownlist <- rbind(WT_DAC_D16_FC_indownlist,WT_BMN_D16_FC_indownlist,WT_DAC_BMN_D16_FC_indownlist)

merge_WT_down_D16 <- merge(Geneset_WT_D16_indownlist,Sequence_WT_D16_indownlist,by.x = "X",by.y = "X")

ggplot(merge_WT_down_D16,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes down-regulated in HIGH-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression change under treatments in WT 786-O cells (Day16)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank())+
  scale_color_manual(values = c("red","blue"))

# KO
Sequence_KO_D16_indownlist <- transform(KO_DAC_D16_FC_indownlist,Genesequence=c(1:657));Sequence_KO_D16_indownlist <- Sequence_KO_D16_indownlist[,c(1,5)]

Geneset_KO_D16_indownlist <- rbind(KO_DAC_D16_FC_indownlist,KO_BMN_D16_FC_indownlist,KO_DAC_BMN_D16_FC_indownlist)

merge_KO_down_D16 <- merge(Geneset_KO_D16_indownlist,Sequence_KO_D16_indownlist,by.x = "X",by.y = "X")

ggplot(merge_KO_down_D16,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes down-regulated in LOW-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression change under treatments in KO 786-O cells (Day16)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank())+
  scale_color_manual(values = c("green","red","blue"))


#Day26
WT_DAC_D26_FC_indownlist <- filtered_edgeRt_WT_D26_DAC_vs_NT;WT_DAC_D26_FC_indownlist <- transform(WT_DAC_D26_FC_indownlist,inpatients=WT_DAC_D26_FC_indownlist$X%in%row.names(SETD2_HIGH_T_vs_N_down_gene));WT_DAC_D26_FC_indownlist <- WT_DAC_D26_FC_indownlist[WT_DAC_D26_FC_indownlist$inpatients==T,];WT_DAC_D26_FC_indownlist <- WT_DAC_D26_FC_indownlist[,c(1,3)];WT_DAC_D26_FC_indownlist <- transform(WT_DAC_D26_FC_indownlist,Treatment="DAC");WT_DAC_D26_FC_indownlist <- transform(WT_DAC_D26_FC_indownlist,Celltype="WT");WT_DAC_D26_FC_indownlist <- WT_DAC_D26_FC_indownlist[order(WT_DAC_D26_FC_indownlist$logFC,decreasing=T),]

WT_BMN_D26_FC_indownlist <- filtered_edgeRt_WT_D26_BMN_vs_NT;WT_BMN_D26_FC_indownlist <- transform(WT_BMN_D26_FC_indownlist,inpatients=WT_BMN_D26_FC_indownlist$X%in%row.names(SETD2_HIGH_T_vs_N_down_gene));WT_BMN_D26_FC_indownlist <- WT_BMN_D26_FC_indownlist[WT_BMN_D26_FC_indownlist$inpatients==T,];WT_BMN_D26_FC_indownlist <- WT_BMN_D26_FC_indownlist[,c(1,3)];WT_BMN_D26_FC_indownlist <- transform(WT_BMN_D26_FC_indownlist,Treatment="BMN");WT_BMN_D26_FC_indownlist <- transform(WT_BMN_D26_FC_indownlist,Celltype="WT")

WT_DAC_BMN_D26_FC_indownlist <- filtered_edgeRt_WT_D26_DAC_BMN_vs_NT;WT_DAC_BMN_D26_FC_indownlist <- transform(WT_DAC_BMN_D26_FC_indownlist,inpatients=WT_DAC_BMN_D26_FC_indownlist$X%in%row.names(SETD2_HIGH_T_vs_N_down_gene));WT_DAC_BMN_D26_FC_indownlist <- WT_DAC_BMN_D26_FC_indownlist[WT_DAC_BMN_D26_FC_indownlist$inpatients==T,];WT_DAC_BMN_D26_FC_indownlist <- WT_DAC_BMN_D26_FC_indownlist[,c(1,3)];WT_DAC_BMN_D26_FC_indownlist <- transform(WT_DAC_BMN_D26_FC_indownlist,Treatment="DAC_BMN");WT_DAC_BMN_D26_FC_indownlist <- transform(WT_DAC_BMN_D26_FC_indownlist,Celltype="WT")

KO_DAC_D26_FC_indownlist <- filtered_edgeRt_KO_D26_DAC_vs_NT;KO_DAC_D26_FC_indownlist <- transform(KO_DAC_D26_FC_indownlist,inpatients=KO_DAC_D26_FC_indownlist$X%in%row.names(SETD2_LOW_T_vs_N_down_gene));KO_DAC_D26_FC_indownlist <- KO_DAC_D26_FC_indownlist[KO_DAC_D26_FC_indownlist$inpatients==T,];KO_DAC_D26_FC_indownlist <- KO_DAC_D26_FC_indownlist[,c(1,3)];KO_DAC_D26_FC_indownlist <- transform(KO_DAC_D26_FC_indownlist,Treatment="DAC");KO_DAC_D26_FC_indownlist <- transform(KO_DAC_D26_FC_indownlist,Celltype="KO");KO_DAC_D26_FC_indownlist <- KO_DAC_D26_FC_indownlist[order(KO_DAC_D26_FC_indownlist$logFC,decreasing=T),]

KO_BMN_D26_FC_indownlist <- filtered_edgeRt_KO_D26_BMN_vs_NT;KO_BMN_D26_FC_indownlist <- transform(KO_BMN_D26_FC_indownlist,inpatients=KO_BMN_D26_FC_indownlist$X%in%row.names(SETD2_LOW_T_vs_N_down_gene));KO_BMN_D26_FC_indownlist <- KO_BMN_D26_FC_indownlist[KO_BMN_D26_FC_indownlist$inpatients==T,];KO_BMN_D26_FC_indownlist <- KO_BMN_D26_FC_indownlist[,c(1,3)];KO_BMN_D26_FC_indownlist <- transform(KO_BMN_D26_FC_indownlist,Treatment="BMN");KO_BMN_D26_FC_indownlist <- transform(KO_BMN_D26_FC_indownlist,Celltype="KO")

KO_DAC_BMN_D26_FC_indownlist <- filtered_edgeRt_KO_D26_DAC_BMN_vs_NT;KO_DAC_BMN_D26_FC_indownlist <- transform(KO_DAC_BMN_D26_FC_indownlist,inpatients=KO_DAC_BMN_D26_FC_indownlist$X%in%row.names(SETD2_LOW_T_vs_N_down_gene));KO_DAC_BMN_D26_FC_indownlist <- KO_DAC_BMN_D26_FC_indownlist[KO_DAC_BMN_D26_FC_indownlist$inpatients==T,];KO_DAC_BMN_D26_FC_indownlist <- KO_DAC_BMN_D26_FC_indownlist[,c(1,3)];KO_DAC_BMN_D26_FC_indownlist <- transform(KO_DAC_BMN_D26_FC_indownlist,Treatment="DAC_BMN");KO_DAC_BMN_D26_FC_indownlist <- transform(KO_DAC_BMN_D26_FC_indownlist,Celltype="KO")

# WT
Sequence_WT_D26_indownlist <- transform(WT_DAC_D26_FC_indownlist,Genesequence=c(1:303));Sequence_WT_D26_indownlist <- Sequence_WT_D26_indownlist[,c(1,5)]

Geneset_WT_D26_indownlist <- rbind(WT_DAC_D26_FC_indownlist,WT_BMN_D26_FC_indownlist,WT_DAC_BMN_D26_FC_indownlist)

merge_WT_down_D26 <- merge(Geneset_WT_D26_indownlist,Sequence_WT_D26_indownlist,by.x = "X",by.y = "X")

ggplot(merge_WT_down_D26,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes down-regulated in HIGH-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression change under treatments in WT 786-O cells (Day26)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank())+
  scale_color_manual(values = c("green","red","blue"))

# KO
Sequence_KO_D26_indownlist <- transform(KO_DAC_D26_FC_indownlist,Genesequence=c(1:653));Sequence_KO_D26_indownlist <- Sequence_KO_D26_indownlist[,c(1,5)]

Geneset_KO_D26_indownlist <- rbind(KO_DAC_D26_FC_indownlist,KO_BMN_D26_FC_indownlist,KO_DAC_BMN_D26_FC_indownlist)

merge_KO_down_D26 <- merge(Geneset_KO_D26_indownlist,Sequence_KO_D26_indownlist,by.x = "X",by.y = "X")

ggplot(merge_KO_down_D26,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes down-regulated in LOW-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression change under treatments in KO 786-O cells (Day26)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank())+
  scale_color_manual(values = c("green","red","blue"))


#make dotplot only contains the reversed genes
#WT
reversed_only_merge_WT_up_D5 <- merge_WT_up_D5[merge_WT_up_D5$logFC<0,]
reversed_only_merge_WT_down_D5 <- merge_WT_down_D5[merge_WT_down_D5$logFC>0,]
reversed_only_merge_WT_up_D16 <- merge_WT_up_D16[merge_WT_up_D16$logFC<0,]
reversed_only_merge_WT_down_D16 <- merge_WT_down_D16[merge_WT_down_D16$logFC>0,]
reversed_only_merge_WT_up_D26 <- merge_WT_up_D26[merge_WT_up_D26$logFC<0,]
reversed_only_merge_WT_down_D26 <- merge_WT_down_D26[merge_WT_down_D26$logFC>0,]

ggplot(reversed_only_merge_WT_up_D5,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes up-regulated in HIGH-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression decreased under treatments in WT 786-O cells (Day5)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank(),plot.title = element_text(size=12))+
  scale_color_manual(values = c("red","blue"))

ggplot(reversed_only_merge_WT_down_D5,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes down-regulated in HIGH-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression increased under treatments in WT 786-O cells (Day5)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank(),plot.title = element_text(size=12))+
  scale_color_manual(values = c("red","blue"))

ggplot(reversed_only_merge_WT_up_D16,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes up-regulated in HIGH-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression decreased under treatments in WT 786-O cells (Day16)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank(),plot.title = element_text(size=12))+
  scale_color_manual(values = c("green","red","blue"))

ggplot(reversed_only_merge_WT_down_D16,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes down-regulated in HIGH-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression increased under treatments in WT 786-O cells (Day16)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank(),plot.title = element_text(size=12))+
  scale_color_manual(values = c("red","blue"))

ggplot(reversed_only_merge_WT_up_D26,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes up-regulated in HIGH-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression decreased under treatments in WT 786-O cells (Day26)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank(),plot.title = element_text(size=12))+
  scale_color_manual(values = c("green","red","blue"))

ggplot(reversed_only_merge_WT_down_D26,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes down-regulated in HIGH-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression increased under treatments in WT 786-O cells (Day26)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank(),plot.title = element_text(size=12))+
  scale_color_manual(values = c("green","red","blue"))

#KO
reversed_only_merge_KO_up_D5 <- merge_KO_up_D5[merge_KO_up_D5$logFC<0,]
reversed_only_merge_KO_down_D5 <- merge_KO_down_D5[merge_KO_down_D5$logFC>0,]
reversed_only_merge_KO_up_D16 <- merge_KO_up_D16[merge_KO_up_D16$logFC<0,]
reversed_only_merge_KO_down_D16 <- merge_KO_down_D16[merge_KO_down_D16$logFC>0,]
reversed_only_merge_KO_up_D26 <- merge_KO_up_D26[merge_KO_up_D26$logFC<0,]
reversed_only_merge_KO_down_D26 <- merge_KO_down_D26[merge_KO_down_D26$logFC>0,]

ggplot(reversed_only_merge_KO_up_D5,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes up-regulated in LOW-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression decreased under treatments in KO 786-O cells (Day5)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank(),plot.title = element_text(size=12))+
  scale_color_manual(values = c("green","red","blue"))

ggplot(reversed_only_merge_KO_down_D5,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes down-regulated in LOW-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression increased under treatments in KO 786-O cells (Day5)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank(),plot.title = element_text(size=12))+
  scale_color_manual(values = c("green","red","blue"))

ggplot(reversed_only_merge_KO_up_D16,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes up-regulated in LOW-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression decreased under treatments in KO 786-O cells (Day16)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank(),plot.title = element_text(size=12))+
  scale_color_manual(values = c("green","red","blue"))

ggplot(reversed_only_merge_KO_down_D16,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes down-regulated in LOW-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression increased under treatments in KO 786-O cells (Day16)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank(),plot.title = element_text(size=12))+
  scale_color_manual(values = c("green","red","blue"))

ggplot(reversed_only_merge_KO_up_D26,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes up-regulated in LOW-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression decreased under treatments in KO 786-O cells (Day26)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank(),plot.title = element_text(size=12))+
  scale_color_manual(values = c("green","red","blue"))

ggplot(reversed_only_merge_KO_down_D26,aes(x=Genesequence,y=logFC,color=Treatment),size=1)+
  xlab("Genes down-regulated in LOW-SETD2 patients (TCGA)")+
  ylab("logFC")+
  labs(title="Gene expression increased under treatments in KO 786-O cells (Day26)",color="Treatment") +
  geom_point(shape=19,alpha = 0.5)+
  facet_grid(Celltype~.)+
  theme_bw()+theme(axis.text.x = element_blank(),plot.title = element_text(size=12))+
  scale_color_manual(values = c("green","red","blue"))

#count the reversed gene numbers of each comparison (|logFC|>1)
#in up list
reversed_number_WT_DAC_D5_FC_inuplist <- WT_DAC_D5_FC_inuplist[WT_DAC_D5_FC_inuplist$logFC<(-1),]#40
reversed_number_WT_DAC_D16_FC_inuplist <- WT_DAC_D16_FC_inuplist[WT_DAC_D16_FC_inuplist$logFC<(-1),]#74
reversed_number_WT_DAC_D26_FC_inuplist <- WT_DAC_D26_FC_inuplist[WT_DAC_D26_FC_inuplist$logFC<(-1),]#59
reversed_number_WT_BMN_D5_FC_inuplist <- WT_BMN_D5_FC_inuplist[WT_BMN_D5_FC_inuplist$logFC<(-1),]#0
reversed_number_WT_BMN_D16_FC_inuplist <- WT_BMN_D16_FC_inuplist[WT_BMN_D16_FC_inuplist$logFC<(-1),]#1
reversed_number_WT_BMN_D26_FC_inuplist <- WT_BMN_D26_FC_inuplist[WT_BMN_D26_FC_inuplist$logFC<(-1),]#39
reversed_number_WT_DAC_BMN_D5_FC_inuplist <- WT_DAC_BMN_D5_FC_inuplist[WT_DAC_BMN_D5_FC_inuplist$logFC<(-1),]#79
reversed_number_WT_DAC_BMN_D16_FC_inuplist <- WT_DAC_BMN_D16_FC_inuplist[WT_DAC_BMN_D16_FC_inuplist$logFC<(-1),]#145
reversed_number_WT_DAC_BMN_D26_FC_inuplist <- WT_DAC_BMN_D26_FC_inuplist[WT_DAC_BMN_D26_FC_inuplist$logFC<(-1),]#119

reversed_number_KO_DAC_D5_FC_inuplist <- KO_DAC_D5_FC_inuplist[KO_DAC_D5_FC_inuplist$logFC<(-1),]#87
reversed_number_KO_DAC_D16_FC_inuplist <- KO_DAC_D16_FC_inuplist[KO_DAC_D16_FC_inuplist$logFC<(-1),]#116
reversed_number_KO_DAC_D26_FC_inuplist <- KO_DAC_D26_FC_inuplist[KO_DAC_D26_FC_inuplist$logFC<(-1),]#165
reversed_number_KO_BMN_D5_FC_inuplist <- KO_BMN_D5_FC_inuplist[KO_BMN_D5_FC_inuplist$logFC<(-1),]#35
reversed_number_KO_BMN_D16_FC_inuplist <- KO_BMN_D16_FC_inuplist[KO_BMN_D16_FC_inuplist$logFC<(-1),]#76
reversed_number_KO_BMN_D26_FC_inuplist <- KO_BMN_D26_FC_inuplist[KO_BMN_D26_FC_inuplist$logFC<(-1),]#202
reversed_number_KO_DAC_BMN_D5_FC_inuplist <- KO_DAC_BMN_D5_FC_inuplist[KO_DAC_BMN_D5_FC_inuplist$logFC<(-1),]#175
reversed_number_KO_DAC_BMN_D16_FC_inuplist <- KO_DAC_BMN_D16_FC_inuplist[KO_DAC_BMN_D16_FC_inuplist$logFC<(-1),]#153
reversed_number_KO_DAC_BMN_D26_FC_inuplist <- KO_DAC_BMN_D26_FC_inuplist[KO_DAC_BMN_D26_FC_inuplist$logFC<(-1),]#108

#in down list
reversed_number_WT_DAC_D5_FC_indownlist <- WT_DAC_D5_FC_indownlist[WT_DAC_D5_FC_indownlist$logFC>1,]#215
reversed_number_WT_DAC_D16_FC_indownlist <- WT_DAC_D16_FC_indownlist[WT_DAC_D16_FC_indownlist$logFC>1,]#220
reversed_number_WT_DAC_D26_FC_indownlist <- WT_DAC_D26_FC_indownlist[WT_DAC_D26_FC_indownlist$logFC>1,]#198
reversed_number_WT_BMN_D5_FC_indownlist <- WT_BMN_D5_FC_indownlist[WT_BMN_D5_FC_indownlist$logFC>1,]#0
reversed_number_WT_BMN_D16_FC_indownlist <- WT_BMN_D16_FC_indownlist[WT_BMN_D16_FC_indownlist$logFC>1,]#0
reversed_number_WT_BMN_D26_FC_indownlist <- WT_BMN_D26_FC_indownlist[WT_BMN_D26_FC_indownlist$logFC>1,]#5
reversed_number_WT_DAC_BMN_D5_FC_indownlist <- WT_DAC_BMN_D5_FC_indownlist[WT_DAC_BMN_D5_FC_indownlist$logFC>1,]#266
reversed_number_WT_DAC_BMN_D16_FC_indownlist <- WT_DAC_BMN_D16_FC_indownlist[WT_DAC_BMN_D16_FC_indownlist$logFC>1,]#217
reversed_number_WT_DAC_BMN_D26_FC_indownlist <- WT_DAC_BMN_D26_FC_indownlist[WT_DAC_BMN_D26_FC_indownlist$logFC>1,]#215

reversed_number_KO_DAC_D5_FC_indownlist <- KO_DAC_D5_FC_indownlist[KO_DAC_D5_FC_indownlist$logFC>1,]#277
reversed_number_KO_DAC_D16_FC_indownlist <- KO_DAC_D16_FC_indownlist[KO_DAC_D16_FC_indownlist$logFC>1,]#279
reversed_number_KO_DAC_D26_FC_indownlist <- KO_DAC_D26_FC_indownlist[KO_DAC_D26_FC_indownlist$logFC>1,]#248
reversed_number_KO_BMN_D5_FC_indownlist <- KO_BMN_D5_FC_indownlist[KO_BMN_D5_FC_indownlist$logFC>1,]#4
reversed_number_KO_BMN_D16_FC_indownlist <- KO_BMN_D16_FC_indownlist[KO_BMN_D16_FC_indownlist$logFC>1,]#6
reversed_number_KO_BMN_D26_FC_indownlist <- KO_BMN_D26_FC_indownlist[KO_BMN_D26_FC_indownlist$logFC>1,]#7
reversed_number_KO_DAC_BMN_D5_FC_indownlist <- KO_DAC_BMN_D5_FC_indownlist[KO_DAC_BMN_D5_FC_indownlist$logFC>1,]#349
reversed_number_KO_DAC_BMN_D16_FC_indownlist <- KO_DAC_BMN_D16_FC_indownlist[KO_DAC_BMN_D16_FC_indownlist$logFC>1,]#231
reversed_number_KO_DAC_BMN_D26_FC_indownlist <- KO_DAC_BMN_D26_FC_indownlist[KO_DAC_BMN_D26_FC_indownlist$logFC>1,]#103

###run pathway analysis
#pathway analysis of reversed_number_WT_DAC_BMN_D5_FC_inuplist
GO_reversed_number_WT_DAC_BMN_D5_FC_inuplist <- reversed_number_WT_DAC_BMN_D5_FC_inuplist$X
GO_geneLists_reversed_number_WT_DAC_BMN_D5_FC_inuplist=data.frame(ensembl_id=GO_reversed_number_WT_DAC_BMN_D5_FC_inuplist)
GO_results_reversed_number_WT_DAC_BMN_D5_FC_inuplist=merge(GO_geneLists_reversed_number_WT_DAC_BMN_D5_FC_inuplist,EG2Ensembl,by.x='ensembl_id',all.x=T)
GO_id_reversed_number_WT_DAC_BMN_D5_FC_inuplist=na.omit(GO_results_reversed_number_WT_DAC_BMN_D5_FC_inuplist$gene_id)
GO_ego_reversed_number_WT_DAC_BMN_D5_FC_inuplist <- enrichGO(OrgDb = "org.Hs.eg.db",gene = GO_id_reversed_number_WT_DAC_BMN_D5_FC_inuplist, ont="BP", pvalueCutoff = 0.05, readable = TRUE)
dotplot(GO_ego_reversed_number_WT_DAC_BMN_D5_FC_inuplist,showCategory=20,title="Enrichment GO reversed_number_WT_DAC_BMN_D5_FC_inuplist")
barplot(GO_ego_reversed_number_WT_DAC_BMN_D5_FC_inuplist,showCategory=20,title="Enrichment GO reversed_number_WT_DAC_BMN_D5_FC_inuplist")

#pathway analysis of reversed_number_KO_DAC_BMN_D5_FC_inuplist
GO_reversed_number_KO_DAC_BMN_D5_FC_inuplist <- reversed_number_KO_DAC_BMN_D5_FC_inuplist$X
GO_geneLists_reversed_number_KO_DAC_BMN_D5_FC_inuplist=data.frame(ensembl_id=GO_reversed_number_KO_DAC_BMN_D5_FC_inuplist)
GO_results_reversed_number_KO_DAC_BMN_D5_FC_inuplist=merge(GO_geneLists_reversed_number_KO_DAC_BMN_D5_FC_inuplist,EG2Ensembl,by.x='ensembl_id',all.x=T)
GO_id_reversed_number_KO_DAC_BMN_D5_FC_inuplist=na.omit(GO_results_reversed_number_KO_DAC_BMN_D5_FC_inuplist$gene_id)
GO_ego_reversed_number_KO_DAC_BMN_D5_FC_inuplist <- enrichGO(OrgDb = "org.Hs.eg.db",gene = GO_id_reversed_number_KO_DAC_BMN_D5_FC_inuplist, ont="BP", pvalueCutoff = 0.05, readable = TRUE)
dotplot(GO_ego_reversed_number_KO_DAC_BMN_D5_FC_inuplist,showCategory=20,title="Enrichment GO reversed_number_KO_DAC_BMN_D5_FC_inuplist")
barplot(GO_ego_reversed_number_KO_DAC_BMN_D5_FC_inuplist,showCategory=20,title="Enrichment GO reversed_number_KO_DAC_BMN_D5_FC_inuplist")

#pathway analysis of reversed_number_WT_DAC_BMN_D5_FC_indownlist
GO_reversed_number_WT_DAC_BMN_D5_FC_indownlist <- reversed_number_WT_DAC_BMN_D5_FC_indownlist$X
GO_geneLists_reversed_number_WT_DAC_BMN_D5_FC_indownlist=data.frame(ensembl_id=GO_reversed_number_WT_DAC_BMN_D5_FC_indownlist)
GO_results_reversed_number_WT_DAC_BMN_D5_FC_indownlist=merge(GO_geneLists_reversed_number_WT_DAC_BMN_D5_FC_indownlist,EG2Ensembl,by.x='ensembl_id',all.x=T)
GO_id_reversed_number_WT_DAC_BMN_D5_FC_indownlist=na.omit(GO_results_reversed_number_WT_DAC_BMN_D5_FC_indownlist$gene_id)
GO_ego_reversed_number_WT_DAC_BMN_D5_FC_indownlist <- enrichGO(OrgDb = "org.Hs.eg.db",gene = GO_id_reversed_number_WT_DAC_BMN_D5_FC_indownlist, ont="BP", pvalueCutoff = 0.05, readable = TRUE)
dotplot(GO_ego_reversed_number_WT_DAC_BMN_D5_FC_indownlist,showCategory=20,title="Enrichment GO reversed_number_WT_DAC_BMN_D5_FC_indownlist")
barplot(GO_ego_reversed_number_WT_DAC_BMN_D5_FC_indownlist,showCategory=20,title="Enrichment GO reversed_number_WT_DAC_BMN_D5_FC_indownlist")

#pathway analysis of reversed_number_KO_DAC_BMN_D5_FC_indownlist
GO_reversed_number_KO_DAC_BMN_D5_FC_indownlist <- reversed_number_KO_DAC_BMN_D5_FC_indownlist$X
GO_geneLists_reversed_number_KO_DAC_BMN_D5_FC_indownlist=data.frame(ensembl_id=GO_reversed_number_KO_DAC_BMN_D5_FC_indownlist)
GO_results_reversed_number_KO_DAC_BMN_D5_FC_indownlist=merge(GO_geneLists_reversed_number_KO_DAC_BMN_D5_FC_indownlist,EG2Ensembl,by.x='ensembl_id',all.x=T)
GO_id_reversed_number_KO_DAC_BMN_D5_FC_indownlist=na.omit(GO_results_reversed_number_KO_DAC_BMN_D5_FC_indownlist$gene_id)
GO_ego_reversed_number_KO_DAC_BMN_D5_FC_indownlist <- enrichGO(OrgDb = "org.Hs.eg.db",gene = GO_id_reversed_number_KO_DAC_BMN_D5_FC_indownlist, ont="BP", pvalueCutoff = 0.05, readable = TRUE)
dotplot(GO_ego_reversed_number_KO_DAC_BMN_D5_FC_indownlist,showCategory=20,title="Enrichment GO reversed_number_KO_DAC_BMN_D5_FC_indownlist")
barplot(GO_ego_reversed_number_KO_DAC_BMN_D5_FC_indownlist,showCategory=20,title="Enrichment GO reversed_number_KO_DAC_BMN_D5_FC_indownlist")
