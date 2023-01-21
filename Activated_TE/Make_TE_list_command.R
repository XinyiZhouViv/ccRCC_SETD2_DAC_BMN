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
