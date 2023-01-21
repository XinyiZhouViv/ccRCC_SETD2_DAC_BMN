#制作TE,SINE,LINE,LTR列表,先做KO vs WT
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




#SINE
SINE_up_KO_vs_WT_DAC <- c(row.names(D5_DAC_KO_vs_WT_up),row.names(D16_DAC_KO_vs_WT_up),row.names(D26_DAC_KO_vs_WT_up));SINE_up_KO_vs_WT_DAC <- transform(SINE_up_KO_vs_WT_DAC,duplication=duplicated(SINE_up_KO_vs_WT_DAC))
SINE_up_KO_vs_WT_DAC <- transform(SINE_up_KO_vs_WT_DAC,isinSINE=SINE_up_KO_vs_WT_DAC$X_data%in%row.names(SINE));SINE_up_KO_vs_WT_DAC <- SINE_up_KO_vs_WT_DAC[SINE_up_KO_vs_WT_DAC$isinSINE==T,]
SINE_up_KO_vs_WT_DAC <- SINE_up_KO_vs_WT_DAC[SINE_up_KO_vs_WT_DAC$duplication==F,];row.names(SINE_up_KO_vs_WT_DAC) <- SINE_up_KO_vs_WT_DAC$X_data
SINE_up_KO_vs_WT_DAC_excludeNT <- transform(SINE_up_KO_vs_WT_DAC,isinNT=row.names(SINE_up_KO_vs_WT_DAC)%in%row.names(NT_KO_vs_WT_up));SINE_up_KO_vs_WT_DAC_excludeNT <- SINE_up_KO_vs_WT_DAC_excludeNT[SINE_up_KO_vs_WT_DAC_excludeNT$isinNT==F,]

SINE_up_KO_vs_WT_BMN <- c(row.names(D5_BMN_KO_vs_WT_up),row.names(D16_BMN_KO_vs_WT_up),row.names(D26_BMN_KO_vs_WT_up));SINE_up_KO_vs_WT_BMN <- transform(SINE_up_KO_vs_WT_BMN,duplication=duplicated(SINE_up_KO_vs_WT_BMN))
SINE_up_KO_vs_WT_BMN <- transform(SINE_up_KO_vs_WT_BMN,isinSINE=SINE_up_KO_vs_WT_BMN$X_data%in%row.names(SINE));SINE_up_KO_vs_WT_BMN <- SINE_up_KO_vs_WT_BMN[SINE_up_KO_vs_WT_BMN$isinSINE==T,]
SINE_up_KO_vs_WT_BMN <- SINE_up_KO_vs_WT_BMN[SINE_up_KO_vs_WT_BMN$duplication==F,];row.names(SINE_up_KO_vs_WT_BMN) <- SINE_up_KO_vs_WT_BMN$X_data
SINE_up_KO_vs_WT_BMN_excludeNT <- transform(SINE_up_KO_vs_WT_BMN,isinNT=row.names(SINE_up_KO_vs_WT_BMN)%in%row.names(NT_KO_vs_WT_up));SINE_up_KO_vs_WT_BMN_excludeNT <- SINE_up_KO_vs_WT_BMN_excludeNT[SINE_up_KO_vs_WT_BMN_excludeNT$isinNT==F,]

SINE_up_KO_vs_WT_DAC_BMN <- c(row.names(D5_DAC_BMN_KO_vs_WT_up),row.names(D16_DAC_BMN_KO_vs_WT_up),row.names(D26_DAC_BMN_KO_vs_WT_up));SINE_up_KO_vs_WT_DAC_BMN <- transform(SINE_up_KO_vs_WT_DAC_BMN,duplication=duplicated(SINE_up_KO_vs_WT_DAC_BMN))
SINE_up_KO_vs_WT_DAC_BMN <- transform(SINE_up_KO_vs_WT_DAC_BMN,isinSINE=SINE_up_KO_vs_WT_DAC_BMN$X_data%in%row.names(SINE));SINE_up_KO_vs_WT_DAC_BMN <- SINE_up_KO_vs_WT_DAC_BMN[SINE_up_KO_vs_WT_DAC_BMN$isinSINE==T,]
SINE_up_KO_vs_WT_DAC_BMN <- SINE_up_KO_vs_WT_DAC_BMN[SINE_up_KO_vs_WT_DAC_BMN$duplication==F,];row.names(SINE_up_KO_vs_WT_DAC_BMN) <- SINE_up_KO_vs_WT_DAC_BMN$X_data
SINE_up_KO_vs_WT_DAC_BMN_excludeNT <- transform(SINE_up_KO_vs_WT_DAC_BMN,isinNT=row.names(SINE_up_KO_vs_WT_DAC_BMN)%in%row.names(NT_KO_vs_WT_up));SINE_up_KO_vs_WT_DAC_BMN_excludeNT <- SINE_up_KO_vs_WT_DAC_BMN_excludeNT[SINE_up_KO_vs_WT_DAC_BMN_excludeNT$isinNT==F,]

SINE_up_KO_vs_WT_ALL <- c(row.names(D5_DAC_KO_vs_WT_up),row.names(D16_DAC_KO_vs_WT_up),row.names(D26_DAC_KO_vs_WT_up),row.names(D5_BMN_KO_vs_WT_up),row.names(D16_BMN_KO_vs_WT_up),row.names(D26_BMN_KO_vs_WT_up),row.names(D5_DAC_BMN_KO_vs_WT_up),row.names(D16_DAC_BMN_KO_vs_WT_up),row.names(D26_DAC_BMN_KO_vs_WT_up));SINE_up_KO_vs_WT_ALL <- transform(SINE_up_KO_vs_WT_ALL,duplication=duplicated(SINE_up_KO_vs_WT_ALL))
SINE_up_KO_vs_WT_ALL <- transform(SINE_up_KO_vs_WT_ALL,isinSINE=SINE_up_KO_vs_WT_ALL$X_data%in%row.names(SINE));SINE_up_KO_vs_WT_ALL <- SINE_up_KO_vs_WT_ALL[SINE_up_KO_vs_WT_ALL$isinSINE==T,]
SINE_up_KO_vs_WT_ALL <- SINE_up_KO_vs_WT_ALL[SINE_up_KO_vs_WT_ALL$duplication==F,];row.names(SINE_up_KO_vs_WT_ALL) <- SINE_up_KO_vs_WT_ALL$X_data
SINE_up_KO_vs_WT_ALL_excludeNT <- transform(SINE_up_KO_vs_WT_ALL,isinNT=row.names(SINE_up_KO_vs_WT_ALL)%in%row.names(NT_KO_vs_WT_up));SINE_up_KO_vs_WT_ALL_excludeNT <- SINE_up_KO_vs_WT_ALL_excludeNT[SINE_up_KO_vs_WT_ALL_excludeNT$isinNT==F,]



#LINE
LINE_up_KO_vs_WT_DAC <- c(row.names(D5_DAC_KO_vs_WT_up),row.names(D16_DAC_KO_vs_WT_up),row.names(D26_DAC_KO_vs_WT_up));LINE_up_KO_vs_WT_DAC <- transform(LINE_up_KO_vs_WT_DAC,duplication=duplicated(LINE_up_KO_vs_WT_DAC))
LINE_up_KO_vs_WT_DAC <- transform(LINE_up_KO_vs_WT_DAC,isinLINE=LINE_up_KO_vs_WT_DAC$X_data%in%row.names(LINE));LINE_up_KO_vs_WT_DAC <- LINE_up_KO_vs_WT_DAC[LINE_up_KO_vs_WT_DAC$isinLINE==T,]
LINE_up_KO_vs_WT_DAC <- LINE_up_KO_vs_WT_DAC[LINE_up_KO_vs_WT_DAC$duplication==F,];row.names(LINE_up_KO_vs_WT_DAC) <- LINE_up_KO_vs_WT_DAC$X_data
LINE_up_KO_vs_WT_DAC_excludeNT <- transform(LINE_up_KO_vs_WT_DAC,isinNT=row.names(LINE_up_KO_vs_WT_DAC)%in%row.names(NT_KO_vs_WT_up));LINE_up_KO_vs_WT_DAC_excludeNT <- LINE_up_KO_vs_WT_DAC_excludeNT[LINE_up_KO_vs_WT_DAC_excludeNT$isinNT==F,]

LINE_up_KO_vs_WT_BMN <- c(row.names(D5_BMN_KO_vs_WT_up),row.names(D16_BMN_KO_vs_WT_up),row.names(D26_BMN_KO_vs_WT_up));LINE_up_KO_vs_WT_BMN <- transform(LINE_up_KO_vs_WT_BMN,duplication=duplicated(LINE_up_KO_vs_WT_BMN))
LINE_up_KO_vs_WT_BMN <- transform(LINE_up_KO_vs_WT_BMN,isinLINE=LINE_up_KO_vs_WT_BMN$X_data%in%row.names(LINE));LINE_up_KO_vs_WT_BMN <- LINE_up_KO_vs_WT_BMN[LINE_up_KO_vs_WT_BMN$isinLINE==T,]
LINE_up_KO_vs_WT_BMN <- LINE_up_KO_vs_WT_BMN[LINE_up_KO_vs_WT_BMN$duplication==F,];row.names(LINE_up_KO_vs_WT_BMN) <- LINE_up_KO_vs_WT_BMN$X_data
LINE_up_KO_vs_WT_BMN_excludeNT <- transform(LINE_up_KO_vs_WT_BMN,isinNT=row.names(LINE_up_KO_vs_WT_BMN)%in%row.names(NT_KO_vs_WT_up));LINE_up_KO_vs_WT_BMN_excludeNT <- LINE_up_KO_vs_WT_BMN_excludeNT[LINE_up_KO_vs_WT_BMN_excludeNT$isinNT==F,]

LINE_up_KO_vs_WT_DAC_BMN <- c(row.names(D5_DAC_BMN_KO_vs_WT_up),row.names(D16_DAC_BMN_KO_vs_WT_up),row.names(D26_DAC_BMN_KO_vs_WT_up));LINE_up_KO_vs_WT_DAC_BMN <- transform(LINE_up_KO_vs_WT_DAC_BMN,duplication=duplicated(LINE_up_KO_vs_WT_DAC_BMN))
LINE_up_KO_vs_WT_DAC_BMN <- transform(LINE_up_KO_vs_WT_DAC_BMN,isinLINE=LINE_up_KO_vs_WT_DAC_BMN$X_data%in%row.names(LINE));LINE_up_KO_vs_WT_DAC_BMN <- LINE_up_KO_vs_WT_DAC_BMN[LINE_up_KO_vs_WT_DAC_BMN$isinLINE==T,]
LINE_up_KO_vs_WT_DAC_BMN <- LINE_up_KO_vs_WT_DAC_BMN[LINE_up_KO_vs_WT_DAC_BMN$duplication==F,];row.names(LINE_up_KO_vs_WT_DAC_BMN) <- LINE_up_KO_vs_WT_DAC_BMN$X_data
LINE_up_KO_vs_WT_DAC_BMN_excludeNT <- transform(LINE_up_KO_vs_WT_DAC_BMN,isinNT=row.names(LINE_up_KO_vs_WT_DAC_BMN)%in%row.names(NT_KO_vs_WT_up));LINE_up_KO_vs_WT_DAC_BMN_excludeNT <- LINE_up_KO_vs_WT_DAC_BMN_excludeNT[LINE_up_KO_vs_WT_DAC_BMN_excludeNT$isinNT==F,]

LINE_up_KO_vs_WT_ALL <- c(row.names(D5_DAC_KO_vs_WT_up),row.names(D16_DAC_KO_vs_WT_up),row.names(D26_DAC_KO_vs_WT_up),row.names(D5_BMN_KO_vs_WT_up),row.names(D16_BMN_KO_vs_WT_up),row.names(D26_BMN_KO_vs_WT_up),row.names(D5_DAC_BMN_KO_vs_WT_up),row.names(D16_DAC_BMN_KO_vs_WT_up),row.names(D26_DAC_BMN_KO_vs_WT_up));LINE_up_KO_vs_WT_ALL <- transform(LINE_up_KO_vs_WT_ALL,duplication=duplicated(LINE_up_KO_vs_WT_ALL))
LINE_up_KO_vs_WT_ALL <- transform(LINE_up_KO_vs_WT_ALL,isinLINE=LINE_up_KO_vs_WT_ALL$X_data%in%row.names(LINE));LINE_up_KO_vs_WT_ALL <- LINE_up_KO_vs_WT_ALL[LINE_up_KO_vs_WT_ALL$isinLINE==T,]
LINE_up_KO_vs_WT_ALL <- LINE_up_KO_vs_WT_ALL[LINE_up_KO_vs_WT_ALL$duplication==F,];row.names(LINE_up_KO_vs_WT_ALL) <- LINE_up_KO_vs_WT_ALL$X_data
LINE_up_KO_vs_WT_ALL_excludeNT <- transform(LINE_up_KO_vs_WT_ALL,isinNT=row.names(LINE_up_KO_vs_WT_ALL)%in%row.names(NT_KO_vs_WT_up));LINE_up_KO_vs_WT_ALL_excludeNT <- LINE_up_KO_vs_WT_ALL_excludeNT[LINE_up_KO_vs_WT_ALL_excludeNT$isinNT==F,]




#LTR
LTR_up_KO_vs_WT_DAC <- c(row.names(D5_DAC_KO_vs_WT_up),row.names(D16_DAC_KO_vs_WT_up),row.names(D26_DAC_KO_vs_WT_up));LTR_up_KO_vs_WT_DAC <- transform(LTR_up_KO_vs_WT_DAC,duplication=duplicated(LTR_up_KO_vs_WT_DAC))
LTR_up_KO_vs_WT_DAC <- transform(LTR_up_KO_vs_WT_DAC,isinLTR=LTR_up_KO_vs_WT_DAC$X_data%in%row.names(LTR));LTR_up_KO_vs_WT_DAC <- LTR_up_KO_vs_WT_DAC[LTR_up_KO_vs_WT_DAC$isinLTR==T,]
LTR_up_KO_vs_WT_DAC <- LTR_up_KO_vs_WT_DAC[LTR_up_KO_vs_WT_DAC$duplication==F,];row.names(LTR_up_KO_vs_WT_DAC) <- LTR_up_KO_vs_WT_DAC$X_data
LTR_up_KO_vs_WT_DAC_excludeNT <- transform(LTR_up_KO_vs_WT_DAC,isinNT=row.names(LTR_up_KO_vs_WT_DAC)%in%row.names(NT_KO_vs_WT_up));LTR_up_KO_vs_WT_DAC_excludeNT <- LTR_up_KO_vs_WT_DAC_excludeNT[LTR_up_KO_vs_WT_DAC_excludeNT$isinNT==F,]

LTR_up_KO_vs_WT_BMN <- c(row.names(D5_BMN_KO_vs_WT_up),row.names(D16_BMN_KO_vs_WT_up),row.names(D26_BMN_KO_vs_WT_up));LTR_up_KO_vs_WT_BMN <- transform(LTR_up_KO_vs_WT_BMN,duplication=duplicated(LTR_up_KO_vs_WT_BMN))
LTR_up_KO_vs_WT_BMN <- transform(LTR_up_KO_vs_WT_BMN,isinLTR=LTR_up_KO_vs_WT_BMN$X_data%in%row.names(LTR));LTR_up_KO_vs_WT_BMN <- LTR_up_KO_vs_WT_BMN[LTR_up_KO_vs_WT_BMN$isinLTR==T,]
LTR_up_KO_vs_WT_BMN <- LTR_up_KO_vs_WT_BMN[LTR_up_KO_vs_WT_BMN$duplication==F,];row.names(LTR_up_KO_vs_WT_BMN) <- LTR_up_KO_vs_WT_BMN$X_data
LTR_up_KO_vs_WT_BMN_excludeNT <- transform(LTR_up_KO_vs_WT_BMN,isinNT=row.names(LTR_up_KO_vs_WT_BMN)%in%row.names(NT_KO_vs_WT_up));LTR_up_KO_vs_WT_BMN_excludeNT <- LTR_up_KO_vs_WT_BMN_excludeNT[LTR_up_KO_vs_WT_BMN_excludeNT$isinNT==F,]

LTR_up_KO_vs_WT_DAC_BMN <- c(row.names(D5_DAC_BMN_KO_vs_WT_up),row.names(D16_DAC_BMN_KO_vs_WT_up),row.names(D26_DAC_BMN_KO_vs_WT_up));LTR_up_KO_vs_WT_DAC_BMN <- transform(LTR_up_KO_vs_WT_DAC_BMN,duplication=duplicated(LTR_up_KO_vs_WT_DAC_BMN))
LTR_up_KO_vs_WT_DAC_BMN <- transform(LTR_up_KO_vs_WT_DAC_BMN,isinLTR=LTR_up_KO_vs_WT_DAC_BMN$X_data%in%row.names(LTR));LTR_up_KO_vs_WT_DAC_BMN <- LTR_up_KO_vs_WT_DAC_BMN[LTR_up_KO_vs_WT_DAC_BMN$isinLTR==T,]
LTR_up_KO_vs_WT_DAC_BMN <- LTR_up_KO_vs_WT_DAC_BMN[LTR_up_KO_vs_WT_DAC_BMN$duplication==F,];row.names(LTR_up_KO_vs_WT_DAC_BMN) <- LTR_up_KO_vs_WT_DAC_BMN$X_data
LTR_up_KO_vs_WT_DAC_BMN_excludeNT <- transform(LTR_up_KO_vs_WT_DAC_BMN,isinNT=row.names(LTR_up_KO_vs_WT_DAC_BMN)%in%row.names(NT_KO_vs_WT_up));LTR_up_KO_vs_WT_DAC_BMN_excludeNT <- LTR_up_KO_vs_WT_DAC_BMN_excludeNT[LTR_up_KO_vs_WT_DAC_BMN_excludeNT$isinNT==F,]

LTR_up_KO_vs_WT_ALL <- c(row.names(D5_DAC_KO_vs_WT_up),row.names(D16_DAC_KO_vs_WT_up),row.names(D26_DAC_KO_vs_WT_up),row.names(D5_BMN_KO_vs_WT_up),row.names(D16_BMN_KO_vs_WT_up),row.names(D26_BMN_KO_vs_WT_up),row.names(D5_DAC_BMN_KO_vs_WT_up),row.names(D16_DAC_BMN_KO_vs_WT_up),row.names(D26_DAC_BMN_KO_vs_WT_up));LTR_up_KO_vs_WT_ALL <- transform(LTR_up_KO_vs_WT_ALL,duplication=duplicated(LTR_up_KO_vs_WT_ALL))
LTR_up_KO_vs_WT_ALL <- transform(LTR_up_KO_vs_WT_ALL,isinLTR=LTR_up_KO_vs_WT_ALL$X_data%in%row.names(LTR));LTR_up_KO_vs_WT_ALL <- LTR_up_KO_vs_WT_ALL[LTR_up_KO_vs_WT_ALL$isinLTR==T,]
LTR_up_KO_vs_WT_ALL <- LTR_up_KO_vs_WT_ALL[LTR_up_KO_vs_WT_ALL$duplication==F,];row.names(LTR_up_KO_vs_WT_ALL) <- LTR_up_KO_vs_WT_ALL$X_data
LTR_up_KO_vs_WT_ALL_excludeNT <- transform(LTR_up_KO_vs_WT_ALL,isinNT=row.names(LTR_up_KO_vs_WT_ALL)%in%row.names(NT_KO_vs_WT_up));LTR_up_KO_vs_WT_ALL_excludeNT <- LTR_up_KO_vs_WT_ALL_excludeNT[LTR_up_KO_vs_WT_ALL_excludeNT$isinNT==F,]

#制作TE,SINE,LINE,LTR列表,再做Treated vs NT
#TE
TE_up_WT_DAC_vs_NT <- c(row.names(WT_D5_DAC_vs_NT_up),row.names(WT_D16_DAC_vs_NT_up),row.names(WT_D26_DAC_vs_NT_up));TE_up_WT_DAC_vs_NT <- transform(TE_up_WT_DAC_vs_NT,duplication=duplicated(TE_up_WT_DAC_vs_NT))
TE_up_WT_DAC_vs_NT <- TE_up_WT_DAC_vs_NT[TE_up_WT_DAC_vs_NT$duplication==F,];row.names(TE_up_WT_DAC_vs_NT) <- TE_up_WT_DAC_vs_NT$X_data
TE_up_KO_DAC_vs_NT <- c(row.names(KO_D5_DAC_vs_NT_up),row.names(KO_D16_DAC_vs_NT_up),row.names(KO_D26_DAC_vs_NT_up));TE_up_KO_DAC_vs_NT <- transform(TE_up_KO_DAC_vs_NT,duplication=duplicated(TE_up_KO_DAC_vs_NT))
TE_up_KO_DAC_vs_NT <- TE_up_KO_DAC_vs_NT[TE_up_KO_DAC_vs_NT$duplication==F,];row.names(TE_up_KO_DAC_vs_NT) <- TE_up_KO_DAC_vs_NT$X_data
TE_up_KO_DAC_vs_NT <- transform(TE_up_KO_DAC_vs_NT,isinWT=row.names(TE_up_KO_DAC_vs_NT)%in%row.names(TE_up_WT_DAC_vs_NT));TE_up_KO_DAC_vs_NT_excludeWT <- TE_up_KO_DAC_vs_NT[TE_up_KO_DAC_vs_NT$isinWT==F,]

TE_up_WT_BMN_vs_NT <- c(row.names(WT_D16_BMN_vs_NT_up),row.names(WT_D26_BMN_vs_NT_up));TE_up_WT_BMN_vs_NT <- transform(TE_up_WT_BMN_vs_NT,duplication=duplicated(TE_up_WT_BMN_vs_NT))
TE_up_WT_BMN_vs_NT <- TE_up_WT_BMN_vs_NT[TE_up_WT_BMN_vs_NT$duplication==F,];row.names(TE_up_WT_BMN_vs_NT) <- TE_up_WT_BMN_vs_NT$X_data
TE_up_KO_BMN_vs_NT <- c(row.names(KO_D5_BMN_vs_NT_up),row.names(KO_D16_BMN_vs_NT_up),row.names(KO_D26_BMN_vs_NT_up));TE_up_KO_BMN_vs_NT <- transform(TE_up_KO_BMN_vs_NT,duplication=duplicated(TE_up_KO_BMN_vs_NT))
TE_up_KO_BMN_vs_NT <- TE_up_KO_BMN_vs_NT[TE_up_KO_BMN_vs_NT$duplication==F,];row.names(TE_up_KO_BMN_vs_NT) <- TE_up_KO_BMN_vs_NT$X_data
TE_up_KO_BMN_vs_NT <- transform(TE_up_KO_BMN_vs_NT,isinWT=row.names(TE_up_KO_BMN_vs_NT)%in%row.names(TE_up_WT_BMN_vs_NT));TE_up_KO_BMN_vs_NT_excludeWT <- TE_up_KO_BMN_vs_NT[TE_up_KO_BMN_vs_NT$isinWT==F,]

TE_up_WT_DAC_BMN_vs_NT <- c(row.names(WT_D5_DAC_BMN_vs_NT_up),row.names(WT_D16_DAC_BMN_vs_NT_up),row.names(WT_D26_DAC_BMN_vs_NT_up));TE_up_WT_DAC_BMN_vs_NT <- transform(TE_up_WT_DAC_BMN_vs_NT,duplication=duplicated(TE_up_WT_DAC_BMN_vs_NT))
TE_up_WT_DAC_BMN_vs_NT <- TE_up_WT_DAC_BMN_vs_NT[TE_up_WT_DAC_BMN_vs_NT$duplication==F,];row.names(TE_up_WT_DAC_BMN_vs_NT) <- TE_up_WT_DAC_BMN_vs_NT$X_data
TE_up_KO_DAC_BMN_vs_NT <- c(row.names(KO_D5_DAC_BMN_vs_NT_up),row.names(KO_D16_DAC_BMN_vs_NT_up),row.names(KO_D26_DAC_BMN_vs_NT_up));TE_up_KO_DAC_BMN_vs_NT <- transform(TE_up_KO_DAC_BMN_vs_NT,duplication=duplicated(TE_up_KO_DAC_BMN_vs_NT))
TE_up_KO_DAC_BMN_vs_NT <- TE_up_KO_DAC_BMN_vs_NT[TE_up_KO_DAC_BMN_vs_NT$duplication==F,];row.names(TE_up_KO_DAC_BMN_vs_NT) <- TE_up_KO_DAC_BMN_vs_NT$X_data
TE_up_KO_DAC_BMN_vs_NT <- transform(TE_up_KO_DAC_BMN_vs_NT,isinWT=row.names(TE_up_KO_DAC_BMN_vs_NT)%in%row.names(TE_up_WT_DAC_BMN_vs_NT));TE_up_KO_DAC_BMN_vs_NT_excludeWT <- TE_up_KO_DAC_BMN_vs_NT[TE_up_KO_DAC_BMN_vs_NT$isinWT==F,]

TE_up_WT_ALL_vs_NT <- c(row.names(TE_up_WT_DAC_vs_NT),row.names(TE_up_WT_BMN_vs_NT),row.names(TE_up_WT_DAC_BMN_vs_NT));TE_up_WT_ALL_vs_NT <- transform(TE_up_WT_ALL_vs_NT,duplication=duplicated(TE_up_WT_ALL_vs_NT))
TE_up_WT_ALL_vs_NT <- TE_up_WT_ALL_vs_NT[TE_up_WT_ALL_vs_NT$duplication==F,];row.names(TE_up_WT_ALL_vs_NT) <- TE_up_WT_ALL_vs_NT$X_data
TE_up_KO_ALL_vs_NT <- c(row.names(TE_up_KO_DAC_vs_NT),row.names(TE_up_KO_BMN_vs_NT),row.names(TE_up_KO_DAC_BMN_vs_NT));TE_up_KO_ALL_vs_NT <- transform(TE_up_KO_ALL_vs_NT,duplication=duplicated(TE_up_KO_ALL_vs_NT))
TE_up_KO_ALL_vs_NT <- TE_up_KO_ALL_vs_NT[TE_up_KO_ALL_vs_NT$duplication==F,];row.names(TE_up_KO_ALL_vs_NT) <- TE_up_KO_ALL_vs_NT$X_data
TE_up_KO_ALL_vs_NT <- transform(TE_up_KO_ALL_vs_NT,isinWT=row.names(TE_up_KO_ALL_vs_NT)%in%row.names(TE_up_WT_ALL_vs_NT));TE_up_KO_ALL_vs_NT_excludeWT <- TE_up_KO_ALL_vs_NT[TE_up_KO_ALL_vs_NT$isinWT==F,]
TE_up_ALL_vs_NT <- c(row.names(TE_up_WT_ALL_vs_NT),row.names(TE_up_KO_ALL_vs_NT));TE_up_ALL_vs_NT <- transform(TE_up_ALL_vs_NT,duplication=duplicated(TE_up_ALL_vs_NT)); TE_up_ALL_vs_NT <- TE_up_ALL_vs_NT[TE_up_ALL_vs_NT$duplication==FALSE,]




#SINE
SINE_up_WT_DAC_vs_NT <- c(row.names(WT_D5_DAC_vs_NT_up),row.names(WT_D16_DAC_vs_NT_up),row.names(WT_D26_DAC_vs_NT_up));SINE_up_WT_DAC_vs_NT <- transform(SINE_up_WT_DAC_vs_NT,duplication=duplicated(SINE_up_WT_DAC_vs_NT))
SINE_up_WT_DAC_vs_NT <- SINE_up_WT_DAC_vs_NT[SINE_up_WT_DAC_vs_NT$duplication==F,];row.names(SINE_up_WT_DAC_vs_NT) <- SINE_up_WT_DAC_vs_NT$X_data
SINE_up_WT_DAC_vs_NT <- transform(SINE_up_WT_DAC_vs_NT,isinSINE=SINE_up_WT_DAC_vs_NT$X_data%in%row.names(SINE));SINE_up_WT_DAC_vs_NT <- SINE_up_WT_DAC_vs_NT[SINE_up_WT_DAC_vs_NT$isinSINE==T,]
SINE_up_KO_DAC_vs_NT <- c(row.names(KO_D5_DAC_vs_NT_up),row.names(KO_D16_DAC_vs_NT_up),row.names(KO_D26_DAC_vs_NT_up));SINE_up_KO_DAC_vs_NT <- transform(SINE_up_KO_DAC_vs_NT,duplication=duplicated(SINE_up_KO_DAC_vs_NT))
SINE_up_KO_DAC_vs_NT <- SINE_up_KO_DAC_vs_NT[SINE_up_KO_DAC_vs_NT$duplication==F,];row.names(SINE_up_KO_DAC_vs_NT) <- SINE_up_KO_DAC_vs_NT$X_data
SINE_up_KO_DAC_vs_NT <- transform(SINE_up_KO_DAC_vs_NT,isinSINE=SINE_up_KO_DAC_vs_NT$X_data%in%row.names(SINE));SINE_up_KO_DAC_vs_NT <- SINE_up_KO_DAC_vs_NT[SINE_up_KO_DAC_vs_NT$isinSINE==T,]
SINE_up_KO_DAC_vs_NT <- transform(SINE_up_KO_DAC_vs_NT,isinWT=row.names(SINE_up_KO_DAC_vs_NT)%in%row.names(SINE_up_WT_DAC_vs_NT));SINE_up_KO_DAC_vs_NT_excludeWT <- SINE_up_KO_DAC_vs_NT[SINE_up_KO_DAC_vs_NT$isinWT==F,]

SINE_up_WT_BMN_vs_NT <- c(row.names(WT_D16_BMN_vs_NT_up),row.names(WT_D26_BMN_vs_NT_up));SINE_up_WT_BMN_vs_NT <- transform(SINE_up_WT_BMN_vs_NT,duplication=duplicated(SINE_up_WT_BMN_vs_NT))
SINE_up_WT_BMN_vs_NT <- SINE_up_WT_BMN_vs_NT[SINE_up_WT_BMN_vs_NT$duplication==F,];row.names(SINE_up_WT_BMN_vs_NT) <- SINE_up_WT_BMN_vs_NT$X_data
SINE_up_WT_BMN_vs_NT <- transform(SINE_up_WT_BMN_vs_NT,isinSINE=SINE_up_WT_BMN_vs_NT$X_data%in%row.names(SINE));SINE_up_WT_BMN_vs_NT <- SINE_up_WT_BMN_vs_NT[SINE_up_WT_BMN_vs_NT$isinSINE==T,]
SINE_up_KO_BMN_vs_NT <- c(row.names(KO_D5_BMN_vs_NT_up),row.names(KO_D16_BMN_vs_NT_up),row.names(KO_D26_BMN_vs_NT_up));SINE_up_KO_BMN_vs_NT <- transform(SINE_up_KO_BMN_vs_NT,duplication=duplicated(SINE_up_KO_BMN_vs_NT))
SINE_up_KO_BMN_vs_NT <- SINE_up_KO_BMN_vs_NT[SINE_up_KO_BMN_vs_NT$duplication==F,];row.names(SINE_up_KO_BMN_vs_NT) <- SINE_up_KO_BMN_vs_NT$X_data
SINE_up_KO_BMN_vs_NT <- transform(SINE_up_KO_BMN_vs_NT,isinSINE=SINE_up_KO_BMN_vs_NT$X_data%in%row.names(SINE));SINE_up_KO_BMN_vs_NT <- SINE_up_KO_BMN_vs_NT[SINE_up_KO_BMN_vs_NT$isinSINE==T,]
SINE_up_KO_BMN_vs_NT <- transform(SINE_up_KO_BMN_vs_NT,isinWT=row.names(SINE_up_KO_BMN_vs_NT)%in%row.names(SINE_up_WT_BMN_vs_NT));SINE_up_KO_BMN_vs_NT_excludeWT <- SINE_up_KO_BMN_vs_NT[SINE_up_KO_BMN_vs_NT$isinWT==F,]

SINE_up_WT_DAC_BMN_vs_NT <- c(row.names(WT_D5_DAC_BMN_vs_NT_up),row.names(WT_D16_DAC_BMN_vs_NT_up),row.names(WT_D26_DAC_BMN_vs_NT_up));SINE_up_WT_DAC_BMN_vs_NT <- transform(SINE_up_WT_DAC_BMN_vs_NT,duplication=duplicated(SINE_up_WT_DAC_BMN_vs_NT))
SINE_up_WT_DAC_BMN_vs_NT <- SINE_up_WT_DAC_BMN_vs_NT[SINE_up_WT_DAC_BMN_vs_NT$duplication==F,];row.names(SINE_up_WT_DAC_BMN_vs_NT) <- SINE_up_WT_DAC_BMN_vs_NT$X_data
SINE_up_WT_DAC_BMN_vs_NT <- transform(SINE_up_WT_DAC_BMN_vs_NT,isinSINE=SINE_up_WT_DAC_BMN_vs_NT$X_data%in%row.names(SINE));SINE_up_WT_DAC_BMN_vs_NT <- SINE_up_WT_DAC_BMN_vs_NT[SINE_up_WT_DAC_BMN_vs_NT$isinSINE==T,]
SINE_up_KO_DAC_BMN_vs_NT <- c(row.names(KO_D5_DAC_BMN_vs_NT_up),row.names(KO_D16_DAC_BMN_vs_NT_up),row.names(KO_D26_DAC_BMN_vs_NT_up));SINE_up_KO_DAC_BMN_vs_NT <- transform(SINE_up_KO_DAC_BMN_vs_NT,duplication=duplicated(SINE_up_KO_DAC_BMN_vs_NT))
SINE_up_KO_DAC_BMN_vs_NT <- SINE_up_KO_DAC_BMN_vs_NT[SINE_up_KO_DAC_BMN_vs_NT$duplication==F,];row.names(SINE_up_KO_DAC_BMN_vs_NT) <- SINE_up_KO_DAC_BMN_vs_NT$X_data
SINE_up_KO_DAC_BMN_vs_NT <- transform(SINE_up_KO_DAC_BMN_vs_NT,isinSINE=SINE_up_KO_DAC_BMN_vs_NT$X_data%in%row.names(SINE));SINE_up_KO_DAC_BMN_vs_NT <- SINE_up_KO_DAC_BMN_vs_NT[SINE_up_KO_DAC_BMN_vs_NT$isinSINE==T,]
SINE_up_KO_DAC_BMN_vs_NT <- transform(SINE_up_KO_DAC_BMN_vs_NT,isinWT=row.names(SINE_up_KO_DAC_BMN_vs_NT)%in%row.names(SINE_up_WT_DAC_BMN_vs_NT));SINE_up_KO_DAC_BMN_vs_NT_excludeWT <- SINE_up_KO_DAC_BMN_vs_NT[SINE_up_KO_DAC_BMN_vs_NT$isinWT==F,]





#LINE
LINE_up_WT_DAC_vs_NT <- c(row.names(WT_D5_DAC_vs_NT_up),row.names(WT_D16_DAC_vs_NT_up),row.names(WT_D26_DAC_vs_NT_up));LINE_up_WT_DAC_vs_NT <- transform(LINE_up_WT_DAC_vs_NT,duplication=duplicated(LINE_up_WT_DAC_vs_NT))
LINE_up_WT_DAC_vs_NT <- LINE_up_WT_DAC_vs_NT[LINE_up_WT_DAC_vs_NT$duplication==F,];row.names(LINE_up_WT_DAC_vs_NT) <- LINE_up_WT_DAC_vs_NT$X_data
LINE_up_WT_DAC_vs_NT <- transform(LINE_up_WT_DAC_vs_NT,isinLINE=LINE_up_WT_DAC_vs_NT$X_data%in%row.names(LINE));LINE_up_WT_DAC_vs_NT <- LINE_up_WT_DAC_vs_NT[LINE_up_WT_DAC_vs_NT$isinLINE==T,]
LINE_up_KO_DAC_vs_NT <- c(row.names(KO_D5_DAC_vs_NT_up),row.names(KO_D16_DAC_vs_NT_up),row.names(KO_D26_DAC_vs_NT_up));LINE_up_KO_DAC_vs_NT <- transform(LINE_up_KO_DAC_vs_NT,duplication=duplicated(LINE_up_KO_DAC_vs_NT))
LINE_up_KO_DAC_vs_NT <- LINE_up_KO_DAC_vs_NT[LINE_up_KO_DAC_vs_NT$duplication==F,];row.names(LINE_up_KO_DAC_vs_NT) <- LINE_up_KO_DAC_vs_NT$X_data
LINE_up_KO_DAC_vs_NT <- transform(LINE_up_KO_DAC_vs_NT,isinLINE=LINE_up_KO_DAC_vs_NT$X_data%in%row.names(LINE));LINE_up_KO_DAC_vs_NT <- LINE_up_KO_DAC_vs_NT[LINE_up_KO_DAC_vs_NT$isinLINE==T,]
LINE_up_KO_DAC_vs_NT <- transform(LINE_up_KO_DAC_vs_NT,isinWT=row.names(LINE_up_KO_DAC_vs_NT)%in%row.names(LINE_up_WT_DAC_vs_NT));LINE_up_KO_DAC_vs_NT_excludeWT <- LINE_up_KO_DAC_vs_NT[LINE_up_KO_DAC_vs_NT$isinWT==F,]

LINE_up_WT_BMN_vs_NT <- c(row.names(WT_D16_BMN_vs_NT_up),row.names(WT_D26_BMN_vs_NT_up));LINE_up_WT_BMN_vs_NT <- transform(LINE_up_WT_BMN_vs_NT,duplication=duplicated(LINE_up_WT_BMN_vs_NT))
LINE_up_WT_BMN_vs_NT <- LINE_up_WT_BMN_vs_NT[LINE_up_WT_BMN_vs_NT$duplication==F,];row.names(LINE_up_WT_BMN_vs_NT) <- LINE_up_WT_BMN_vs_NT$X_data
LINE_up_WT_BMN_vs_NT <- transform(LINE_up_WT_BMN_vs_NT,isinLINE=LINE_up_WT_BMN_vs_NT$X_data%in%row.names(LINE));LINE_up_WT_BMN_vs_NT <- LINE_up_WT_BMN_vs_NT[LINE_up_WT_BMN_vs_NT$isinLINE==T,]
LINE_up_KO_BMN_vs_NT <- c(row.names(KO_D5_BMN_vs_NT_up),row.names(KO_D16_BMN_vs_NT_up),row.names(KO_D26_BMN_vs_NT_up));LINE_up_KO_BMN_vs_NT <- transform(LINE_up_KO_BMN_vs_NT,duplication=duplicated(LINE_up_KO_BMN_vs_NT))
LINE_up_KO_BMN_vs_NT <- LINE_up_KO_BMN_vs_NT[LINE_up_KO_BMN_vs_NT$duplication==F,];row.names(LINE_up_KO_BMN_vs_NT) <- LINE_up_KO_BMN_vs_NT$X_data
LINE_up_KO_BMN_vs_NT <- transform(LINE_up_KO_BMN_vs_NT,isinLINE=LINE_up_KO_BMN_vs_NT$X_data%in%row.names(LINE));LINE_up_KO_BMN_vs_NT <- LINE_up_KO_BMN_vs_NT[LINE_up_KO_BMN_vs_NT$isinLINE==T,]
LINE_up_KO_BMN_vs_NT <- transform(LINE_up_KO_BMN_vs_NT,isinWT=row.names(LINE_up_KO_BMN_vs_NT)%in%row.names(LINE_up_WT_BMN_vs_NT));LINE_up_KO_BMN_vs_NT_excludeWT <- LINE_up_KO_BMN_vs_NT[LINE_up_KO_BMN_vs_NT$isinWT==F,]

LINE_up_WT_DAC_BMN_vs_NT <- c(row.names(WT_D5_DAC_BMN_vs_NT_up),row.names(WT_D16_DAC_BMN_vs_NT_up),row.names(WT_D26_DAC_BMN_vs_NT_up));LINE_up_WT_DAC_BMN_vs_NT <- transform(LINE_up_WT_DAC_BMN_vs_NT,duplication=duplicated(LINE_up_WT_DAC_BMN_vs_NT))
LINE_up_WT_DAC_BMN_vs_NT <- LINE_up_WT_DAC_BMN_vs_NT[LINE_up_WT_DAC_BMN_vs_NT$duplication==F,];row.names(LINE_up_WT_DAC_BMN_vs_NT) <- LINE_up_WT_DAC_BMN_vs_NT$X_data
LINE_up_WT_DAC_BMN_vs_NT <- transform(LINE_up_WT_DAC_BMN_vs_NT,isinLINE=LINE_up_WT_DAC_BMN_vs_NT$X_data%in%row.names(LINE));LINE_up_WT_DAC_BMN_vs_NT <- LINE_up_WT_DAC_BMN_vs_NT[LINE_up_WT_DAC_BMN_vs_NT$isinLINE==T,]
LINE_up_KO_DAC_BMN_vs_NT <- c(row.names(KO_D5_DAC_BMN_vs_NT_up),row.names(KO_D16_DAC_BMN_vs_NT_up),row.names(KO_D26_DAC_BMN_vs_NT_up));LINE_up_KO_DAC_BMN_vs_NT <- transform(LINE_up_KO_DAC_BMN_vs_NT,duplication=duplicated(LINE_up_KO_DAC_BMN_vs_NT))
LINE_up_KO_DAC_BMN_vs_NT <- LINE_up_KO_DAC_BMN_vs_NT[LINE_up_KO_DAC_BMN_vs_NT$duplication==F,];row.names(LINE_up_KO_DAC_BMN_vs_NT) <- LINE_up_KO_DAC_BMN_vs_NT$X_data
LINE_up_KO_DAC_BMN_vs_NT <- transform(LINE_up_KO_DAC_BMN_vs_NT,isinLINE=LINE_up_KO_DAC_BMN_vs_NT$X_data%in%row.names(LINE));LINE_up_KO_DAC_BMN_vs_NT <- LINE_up_KO_DAC_BMN_vs_NT[LINE_up_KO_DAC_BMN_vs_NT$isinLINE==T,]
LINE_up_KO_DAC_BMN_vs_NT <- transform(LINE_up_KO_DAC_BMN_vs_NT,isinWT=row.names(LINE_up_KO_DAC_BMN_vs_NT)%in%row.names(LINE_up_WT_DAC_BMN_vs_NT));LINE_up_KO_DAC_BMN_vs_NT_excludeWT <- LINE_up_KO_DAC_BMN_vs_NT[LINE_up_KO_DAC_BMN_vs_NT$isinWT==F,]





#LTR
LTR_up_WT_DAC_vs_NT <- c(row.names(WT_D5_DAC_vs_NT_up),row.names(WT_D16_DAC_vs_NT_up),row.names(WT_D26_DAC_vs_NT_up));LTR_up_WT_DAC_vs_NT <- transform(LTR_up_WT_DAC_vs_NT,duplication=duplicated(LTR_up_WT_DAC_vs_NT))
LTR_up_WT_DAC_vs_NT <- LTR_up_WT_DAC_vs_NT[LTR_up_WT_DAC_vs_NT$duplication==F,];row.names(LTR_up_WT_DAC_vs_NT) <- LTR_up_WT_DAC_vs_NT$X_data
LTR_up_WT_DAC_vs_NT <- transform(LTR_up_WT_DAC_vs_NT,isinLTR=LTR_up_WT_DAC_vs_NT$X_data%in%row.names(LTR));LTR_up_WT_DAC_vs_NT <- LTR_up_WT_DAC_vs_NT[LTR_up_WT_DAC_vs_NT$isinLTR==T,]
LTR_up_KO_DAC_vs_NT <- c(row.names(KO_D5_DAC_vs_NT_up),row.names(KO_D16_DAC_vs_NT_up),row.names(KO_D26_DAC_vs_NT_up));LTR_up_KO_DAC_vs_NT <- transform(LTR_up_KO_DAC_vs_NT,duplication=duplicated(LTR_up_KO_DAC_vs_NT))
LTR_up_KO_DAC_vs_NT <- LTR_up_KO_DAC_vs_NT[LTR_up_KO_DAC_vs_NT$duplication==F,];row.names(LTR_up_KO_DAC_vs_NT) <- LTR_up_KO_DAC_vs_NT$X_data
LTR_up_KO_DAC_vs_NT <- transform(LTR_up_KO_DAC_vs_NT,isinLTR=LTR_up_KO_DAC_vs_NT$X_data%in%row.names(LTR));LTR_up_KO_DAC_vs_NT <- LTR_up_KO_DAC_vs_NT[LTR_up_KO_DAC_vs_NT$isinLTR==T,]
LTR_up_KO_DAC_vs_NT <- transform(LTR_up_KO_DAC_vs_NT,isinWT=row.names(LTR_up_KO_DAC_vs_NT)%in%row.names(LTR_up_WT_DAC_vs_NT));LTR_up_KO_DAC_vs_NT_excludeWT <- LTR_up_KO_DAC_vs_NT[LTR_up_KO_DAC_vs_NT$isinWT==F,]

LTR_up_WT_BMN_vs_NT <- c(row.names(WT_D16_BMN_vs_NT_up),row.names(WT_D26_BMN_vs_NT_up));LTR_up_WT_BMN_vs_NT <- transform(LTR_up_WT_BMN_vs_NT,duplication=duplicated(LTR_up_WT_BMN_vs_NT))
LTR_up_WT_BMN_vs_NT <- LTR_up_WT_BMN_vs_NT[LTR_up_WT_BMN_vs_NT$duplication==F,];row.names(LTR_up_WT_BMN_vs_NT) <- LTR_up_WT_BMN_vs_NT$X_data
LTR_up_WT_BMN_vs_NT <- transform(LTR_up_WT_BMN_vs_NT,isinLTR=LTR_up_WT_BMN_vs_NT$X_data%in%row.names(LTR));LTR_up_WT_BMN_vs_NT <- LTR_up_WT_BMN_vs_NT[LTR_up_WT_BMN_vs_NT$isinLTR==T,]
LTR_up_KO_BMN_vs_NT <- c(row.names(KO_D5_BMN_vs_NT_up),row.names(KO_D16_BMN_vs_NT_up),row.names(KO_D26_BMN_vs_NT_up));LTR_up_KO_BMN_vs_NT <- transform(LTR_up_KO_BMN_vs_NT,duplication=duplicated(LTR_up_KO_BMN_vs_NT))
LTR_up_KO_BMN_vs_NT <- LTR_up_KO_BMN_vs_NT[LTR_up_KO_BMN_vs_NT$duplication==F,];row.names(LTR_up_KO_BMN_vs_NT) <- LTR_up_KO_BMN_vs_NT$X_data
LTR_up_KO_BMN_vs_NT <- transform(LTR_up_KO_BMN_vs_NT,isinLTR=LTR_up_KO_BMN_vs_NT$X_data%in%row.names(LTR));LTR_up_KO_BMN_vs_NT <- LTR_up_KO_BMN_vs_NT[LTR_up_KO_BMN_vs_NT$isinLTR==T,]
LTR_up_KO_BMN_vs_NT <- transform(LTR_up_KO_BMN_vs_NT,isinWT=row.names(LTR_up_KO_BMN_vs_NT)%in%row.names(LTR_up_WT_BMN_vs_NT));LTR_up_KO_BMN_vs_NT_excludeWT <- LTR_up_KO_BMN_vs_NT[LTR_up_KO_BMN_vs_NT$isinWT==F,]

LTR_up_WT_DAC_BMN_vs_NT <- c(row.names(WT_D5_DAC_BMN_vs_NT_up),row.names(WT_D16_DAC_BMN_vs_NT_up),row.names(WT_D26_DAC_BMN_vs_NT_up));LTR_up_WT_DAC_BMN_vs_NT <- transform(LTR_up_WT_DAC_BMN_vs_NT,duplication=duplicated(LTR_up_WT_DAC_BMN_vs_NT))
LTR_up_WT_DAC_BMN_vs_NT <- LTR_up_WT_DAC_BMN_vs_NT[LTR_up_WT_DAC_BMN_vs_NT$duplication==F,];row.names(LTR_up_WT_DAC_BMN_vs_NT) <- LTR_up_WT_DAC_BMN_vs_NT$X_data
LTR_up_WT_DAC_BMN_vs_NT <- transform(LTR_up_WT_DAC_BMN_vs_NT,isinLTR=LTR_up_WT_DAC_BMN_vs_NT$X_data%in%row.names(LTR));LTR_up_WT_DAC_BMN_vs_NT <- LTR_up_WT_DAC_BMN_vs_NT[LTR_up_WT_DAC_BMN_vs_NT$isinLTR==T,]
LTR_up_KO_DAC_BMN_vs_NT <- c(row.names(KO_D5_DAC_BMN_vs_NT_up),row.names(KO_D16_DAC_BMN_vs_NT_up),row.names(KO_D26_DAC_BMN_vs_NT_up));LTR_up_KO_DAC_BMN_vs_NT <- transform(LTR_up_KO_DAC_BMN_vs_NT,duplication=duplicated(LTR_up_KO_DAC_BMN_vs_NT))
LTR_up_KO_DAC_BMN_vs_NT <- LTR_up_KO_DAC_BMN_vs_NT[LTR_up_KO_DAC_BMN_vs_NT$duplication==F,];row.names(LTR_up_KO_DAC_BMN_vs_NT) <- LTR_up_KO_DAC_BMN_vs_NT$X_data
LTR_up_KO_DAC_BMN_vs_NT <- transform(LTR_up_KO_DAC_BMN_vs_NT,isinLTR=LTR_up_KO_DAC_BMN_vs_NT$X_data%in%row.names(LTR));LTR_up_KO_DAC_BMN_vs_NT <- LTR_up_KO_DAC_BMN_vs_NT[LTR_up_KO_DAC_BMN_vs_NT$isinLTR==T,]
LTR_up_KO_DAC_BMN_vs_NT <- transform(LTR_up_KO_DAC_BMN_vs_NT,isinWT=row.names(LTR_up_KO_DAC_BMN_vs_NT)%in%row.names(LTR_up_WT_DAC_BMN_vs_NT));LTR_up_KO_DAC_BMN_vs_NT_excludeWT <- LTR_up_KO_DAC_BMN_vs_NT[LTR_up_KO_DAC_BMN_vs_NT$isinWT==F,]
