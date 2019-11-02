#' ---
#' title: "Differential abundance testing for Hu et. al."
#' author: "F.M. Kerckhof"
#' date: "`r format(Sys.time(), '%F')`"
#' ---



#+ setup, include = FALSE
#### load required packages ####
library(readxl)
library(dplyr)
library(DESeq2)
library(knitr)
library(ggplot2)
library(gtools)
library(scales)
library(ggthemes)
#font_import()
loadfonts(device = "win")
#### User defined functions ####
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#source("Highqualgraphr.R")
source("DataLoader.R")
absconc <- read.csv("Source_Data/Xiaona_absolute_counts.csv")
absconc.relab <- read.csv("Source_Data/Xiaona_absolute_counts.csv")


#+ spacer

#' # Loading & preformatting the data
#' First, all data was loaded and preformated/filtered in different ways for 
#' comparative purposes.

#+ dataloader, fig.width=6, fig.height=7, cache=T
#### load the data ####
meta <- metadata.no.inoculum
shared <- shared.no.inoculum

#### preformat the data ###
sh.x <- shared
# rownames(sh.x) <- shared$OTU
#remove absoute singletons
sh.x.ns <- sh.x[which(rowSums(sh.x)!=1),]
meta <- data.frame(meta,Factor3=factor(paste(meta$Factor1,meta$Factor2,sep="_"))) # to have a one-factor DESeq design
#### prevalence filtering ####
## prevalence = (fraction of samples in which an OTU is observed minimum 1 time)
minobs=1
prevalence <- apply(as.matrix(sh.x.ns),1,function(x,minobs){sum(x>=minobs)},minobs)/ncol(sh.x.ns)
prevalencefilter <- prevalence>0.05
sh.x.prevfilt <- sh.x.ns[prevalencefilter,]
##Read counts should exceed 0.5 times the number of samples
sh.x.prevfilt <- sh.x.prevfilt[rowSums(sh.x.prevfilt)>0.5*ncol(sh.x.prevfilt),]
sh.x.prevfilt.prop <- sweep(sh.x.prevfilt,2,colSums(sh.x.prevfilt),'/')
# deseq normalise ==> very dependent upon design

#+ spacer

#' # Estimating overdispersion and variance stabilisation
#' Next, the prevalence-filtered data was normalized and variance-stabilized.


#+ deseqformatting
deseqdata <- as.matrix(sh.x.prevfilt +1)
deseqdata <- DESeqDataSetFromMatrix(deseqdata,colData=meta,design= ~ Factor3)
deseqdata <- estimateSizeFactors(deseqdata)
sharedfiltereddeseq <- counts(deseqdata,normalized=TRUE)
deseqdata <- estimateDispersions(deseqdata,fitType="local")
sharedfiltereddeseqvartransf <- varianceStabilizingTransformation(deseqdata, blind = FALSE)
sharedfiltereddeseqvartransfmatrix <- assay(sharedfiltereddeseqvartransf)
wnwndeseqvartransf <- getVarianceStabilizedData(deseqdata)
#jaccard_deseq <- vegdist(t(sharedfiltereddeseq),method="jaccard",binary="FALSE") 


#+ spacer3

#' # Differential abundance testing
#' Ultimately, the transformed data was used to test for differential abundance.
#' Multiple testing correction controls the false discovery rate (FDR) using the
#' method proposed by Benjamini & Yekutieli (2001)


#+ DEtestdeseq, cache=TRUE
testdiff <- DESeq(deseqdata,test = "Wald",fitType = "local")
testdiff.res <- results(testdiff,pAdjustMethod = "BY") 
testdiff.res.ordered <- testdiff.res[order(testdiff.res$padj),]
summary(testdiff.res)
signotu <- testdiff.res$padj<0.05
otunames <- row.names(sh.x.prevfilt)[signotu==TRUE]
resdf <- data.frame(sh.x.prevfilt,sigdif=signotu)
write.csv2(resdf,"differentialabundancexiaona.csv")
# with comments from ?results
testdiff.res5 <- results(testdiff,alpha=0.05,pAdjustMethod = "BY")
testdiff.res.ordered5 <- testdiff.res5[order(testdiff.res5$padj),]
summary(testdiff.res5)
signotu5 <- testdiff.res5$padj<0.05
resdf5 <- data.frame(sh.x.prevfilt,sigdif=signotu5,padjusted=testdiff.res5$padj)
write.csv(resdf5,"differentialxiaona5.csv")
#+ spacer 4

#' Hence it appears that `r length(otunames)` OTUs are significantly different among the different sources and nitrogen conditions among
#' the total of `r nrow(sh.x.prevfilt)` OTUs. A csv was created with this information included.


#+ spacer 4

#' # Visualisation all OTU DE
#' Only significant OTUs are plotted below
#' 

pwvals <- results(testdiff,contrast = c("Factor3","BS_N2","BS_NH4"),pAdjustMethod = "BY")$padj
ggcols2 <- gg_color_hue(2)

saveotu.sig <- data.frame()
significantotus <- rownames(resdf5[!is.na(resdf5$sigdif)&resdf5$sigdif==TRUE,])
for(i in significantotus){
  diffotu <- plotCounts(testdiff,gene=i,intgroup = "Factor3",returnData = TRUE,transform=FALSE)
  diffotu <- cbind(diffotu,OTU=i)
  saveotu.sig <- rbind(saveotu.sig,diffotu)
}


rankedsig <- saveotu.sig %>%  group_by(OTU) %>% summarise(med=median(count)) %>% arrange(desc(med)) %>% select(OTU)


#### generate separate ggplot2 objects to be able to extract outliers ####
plotsig.vals <- ggplot(aes(y=count,x=Factor3,color=Factor3,shape=Factor3),
                       data=saveotu.sig) + 
  geom_boxplot(alpha=0.5,outlier.shape = 8,outlier.size = 3,
               outlier.alpha = 1) +
  geom_jitter(position=position_jitter(width=.3,
                                       height=0),
              size=4,alpha=1) +
  scale_y_continuous(labels=scientific,trans='log10',
                     breaks = trans_breaks('log10',
                                           function(x) 10^x)) +
  ylab("Normalised log transformed counts")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size=18)) +
  scale_x_discrete(breaks=c("BS_N2","HE_N2","SS_N2",
                            "BS_NH4","HE_NH4","SS_NH4"),
                   labels=c("BS_N2","HE_N2","SS_N2",
                            "BS_NH4","HE_NH4","SS_NH4")) +
  facet_grid(~OTU,scales="free")

saveotu.sig$Factor3<- factor(saveotu.sig$Factor3,
                             levels=c("BS_N2","HE_N2","SS_N2",
                                      "BS_NH4","HE_NH4","SS_NH4"))

saveotu.sig.XH <- saveotu.sig %>% 
  dplyr::filter(OTU %in% c("Otu00001","Otu00002","Otu00003"))

saveotu.sig.XH$OTUclass <- factor(saveotu.sig.XH$OTU)

levels(saveotu.sig.XH$OTUclass)[levels(saveotu.sig.XH$OTUclass)=="Otu00001"] <-
  "Otu00001 (Azonexus)"

levels(saveotu.sig.XH$OTUclass)[levels(saveotu.sig.XH$OTUclass)=="Otu00002"] <-
  "Otu00002 (Comamonadaceae)"

levels(saveotu.sig.XH$OTUclass)[levels(saveotu.sig.XH$OTUclass)=="Otu00003"] <-
  "Otu00003 (Xanthobacter)"

plotsig.xh2 <- ggplot(aes(y=count,x=Factor3,color=Factor3,shape=Factor3),
                      data=saveotu.sig.XH) + 
  geom_boxplot(alpha=0.5,outlier.shape = 8,outlier.size = 3,
               outlier.alpha = 1) +
  geom_jitter(position=position_jitter(width=.3,
                                       height=0),
              size=2,alpha=1) +
  scale_y_continuous(labels=scientific,trans='log10',
                     breaks = trans_breaks('log10',
                                           function(x) 10^x)) +
  ylab("Normalised log transformed counts")+
  scale_color_discrete(name="Treatment",
                       breaks=c("BS_N2","HE_N2","SS_N2",
                                "BS_NH4","HE_NH4","SS_NH4"),
                       labels=expression(BS-N[2],HE-N[2],SS-N[2],
                                         BS-NH[4]^"+",HE-NH[4]^"+",
                                         SS-NH[4]^"+")) +
  scale_shape_discrete(name="Treatment",
                       breaks=c("BS_N2","HE_N2","SS_N2",
                                "BS_NH4","HE_NH4","SS_NH4"),
                       labels=expression(BS-N[2],HE-N[2],SS-N[2],
                                         BS-NH[4]^"+",HE-NH[4]^"+",
                                         SS-NH[4]^"+")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size=10),
        legend.title= element_blank(),
        legend.text.align = 0,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_grid(~OTUclass,scales="free")

ggsave(filename="XiaonaDEsell.tiff",device="tiff",width=7,
       height=4,units = "in",dpi = 300)

### Boxplot for the same with the counts from the data 
absco.long <- absconc %>% group_by(X) %>% 
  tidyr::gather("Sample_name","Count",-X) %>% 
  filter(X %in% c("Otu00001","Otu00002","Otu00003"))

absco.long$Sample_type <- "BS N2"
absco.long$Sample_type[grep("HE[0-9].*",absco.long$Sample_name)]<-"HE N2"
absco.long$Sample_type[grep("SS[0-9].*",absco.long$Sample_name)]<-"SS N2"
absco.long$Sample_type[grep("BSN[0-9].*",absco.long$Sample_name)]<-"BS NH4"
absco.long$Sample_type[grep("HEN[0-9].*",absco.long$Sample_name)]<-"HE NH4"
absco.long$Sample_type[grep("SSN[0-9].*",absco.long$Sample_name)]<-"SS NH4"

absco.long$Sample_type <- factor(absco.long$Sample_type,
                                 levels=c("BS N2","HE N2","SS N2",
                                          "BS NH4","HE NH4","SS NH4"))
absco.long$OTUclass <- factor(absco.long$X)

levels(absco.long$OTUclass)[levels(absco.long$OTUclass)=="Otu00001"] <-
  "Otu00001 (Azonexus)"

levels(absco.long$OTUclass)[levels(absco.long$OTUclass)=="Otu00002"] <-
  "Otu00002 (Comamonadaceae)"

levels(absco.long$OTUclass)[levels(absco.long$OTUclass)=="Otu00003"] <-
  "Otu00003 (Xanthobacter)"

absco.plot <- ggplot(aes(x=Sample_type,y=Count,shape=Sample_type,
                         color=Sample_type),
                     data=absco.long)+ 
  geom_boxplot(alpha=0.5,outlier.shape = 8,outlier.size = 3,
               outlier.alpha = 1) +
  geom_jitter(position=position_jitter(width=.3,
                                       height=0),
              size=2,alpha=1) +
  ylab("Absolute taxon abundance per mL")+
  scale_color_discrete(name="Treatment",
                       breaks=c("BS N2","HE N2","SS N2",
                                "BS NH4","HE NH4","SS NH4"),
                       labels=expression(BS-N[2],HE-N[2],SS-N[2],
                                         BS-NH[4]^"+",HE-NH[4]^"+",
                                         SS-NH[4]^"+")) +
  scale_shape_discrete(name="Treatment",
                       breaks=c("BS N2","HE N2","SS N2",
                                "BS NH4","HE NH4","SS NH4"),
                       labels=expression(BS-N[2],HE-N[2],SS-N[2],
                                         BS-NH[4]^"+",HE-NH[4]^"+",
                                         SS-NH[4]^"+")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size=10),
        legend.title= element_blank(),
        legend.text.align = 0,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_grid(~OTUclass,scales="free")

ggsave(filename="Xiaonaabscosell.tiff",device="tiff",width=7,
       height=4,units = "in",dpi = 300)

pdf("xiaonaDE.pdf")
for(i in significantotus){
  otusel <- saveotu.sig[saveotu.sig$OTU==i,]
  plotsig.vals <- ggplot(aes(y=count,x=Factor3,color=Factor3,shape=Factor3),
                         data=otusel) + 
    geom_boxplot(alpha=0.5,outlier.shape = 8,outlier.size = 3,
                 outlier.alpha = 1) +
    geom_jitter(position=position_jitter(width=.3,
                                         height=0),
                size=4,alpha=1) +
    scale_y_continuous(labels=scientific,trans='log10',
                       breaks = trans_breaks('log10',
                                             function(x) 10^x)) +
    ylab("Normalised log transformed counts")+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          text = element_text(size=18)) +
    scale_x_discrete(breaks=c("BS_N2","HE_N2","SS_N2",
                              "BS_NH4","HE_NH4","SS_NH4"),
                     labels=c("BS_N2","HE_N2","SS_N2",
                              "BS_NH4","HE_NH4","SS_NH4")) +
    facet_grid(~OTU,scales="free")  
  print(plotsig.vals)
  
}
dev.off()




#### do DE testing in subsampled datasets for each phase ####
#### 1. Phase 1 ####
shared.phase1 <- shared.t.ns %>% 
  dplyr::select(which(grepl(pattern = "[NES]1",colnames(shared.t.ns)))) %>%  
  dplyr::select(-contains("0"))
metadata.phase1 <- metadata[metadata$SampleName %in% colnames(shared.phase1),]

meta <- metadata.phase1
shared <- shared.phase1
sh.x <- shared
sh.x.ns <- sh.x[which(rowSums(sh.x)!=1),]
meta <- data.frame(meta,Factor3=factor(paste(meta$Factor1,meta$Factor2))) # to have a one-factor DESeq design
## prevalence = (fraction of samples in which an OTU is observed minimum 1 time)
minobs=1
prevalence <- apply(as.matrix(sh.x.ns),1,function(x,minobs){sum(x>=minobs)},minobs)/ncol(sh.x.ns)
prevalencefilter <- prevalence>0.05
sh.x.prevfilt <- sh.x.ns[prevalencefilter,]
##Read counts should exceed 0.5 times the number of samples
sh.x.prevfilt <- sh.x.prevfilt[rowSums(sh.x.prevfilt)>0.5*ncol(sh.x.prevfilt),]
sh.x.prevfilt.prop <- sweep(sh.x.prevfilt,2,colSums(sh.x.prevfilt),'/')
deseqdata <- as.matrix(sh.x.prevfilt +1)
deseqdata <- DESeqDataSetFromMatrix(deseqdata,colData=meta,design= ~ Factor3)
deseqdata <- estimateSizeFactors(deseqdata)
sharedfiltereddeseq <- counts(deseqdata,normalized=TRUE)
deseqdata <- estimateDispersions(deseqdata,fitType="local")
sharedfiltereddeseqvartransf <- varianceStabilizingTransformation(deseqdata, blind = FALSE)
sharedfiltereddeseqvartransfmatrix <- assay(sharedfiltereddeseqvartransf)
wnwndeseqvartransf <- getVarianceStabilizedData(deseqdata)
testdiff <- DESeq(deseqdata,test = "Wald",fitType = "local")
testdiff.res5 <- results(testdiff,alpha=0.05,pAdjustMethod = "BY")
testdiff.res.ordered5 <- testdiff.res5[order(testdiff.res5$padj),]
summary(testdiff.res5)
signotu5 <- testdiff.res5$padj<0.05
resdf5 <- data.frame(sh.x.prevfilt,sigdif=signotu5,padjusted=testdiff.res5$padj)
write.csv(resdf5,"differentialxiaona5_phase1.csv")
pwvals <- results(testdiff,contrast = c("Factor3","BS Nfree","BS Nitrogen"),pAdjustMethod = "BY")$padj
saveotu.sig <- data.frame()
significantotus <- rownames(resdf5[!is.na(resdf5$sigdif)&resdf5$sigdif==TRUE,])
for(i in significantotus){
  diffotu <- plotCounts(testdiff,gene=i,intgroup = "Factor3",returnData = TRUE,transform=FALSE)
  diffotu <- cbind(diffotu,OTU=i)
  saveotu.sig <- rbind(saveotu.sig,diffotu)
}
saveotu.sig <- data.frame(saveotu.sig,Source=factor(substr(as.character(saveotu.sig$Factor3),0,2)))
rankedsig <- saveotu.sig %>%  group_by(OTU) %>% summarise(med=median(count)) %>% arrange(desc(med)) %>% select(OTU)
plotsig.vals <- ggplot(aes(y=count,x=Factor3,color=Factor3,shape=Factor3),
                       data=saveotu.sig) + 
  geom_boxplot(alpha=0.5,outlier.shape = 8,outlier.size = 3,
               outlier.alpha = 1) +
  geom_jitter(position=position_jitter(width=.3,
                                       height=0),
              size=4,alpha=1) +
  scale_y_continuous(labels=scientific,trans='log10',
                     breaks = trans_breaks('log10',
                                           function(x) 10^x)) +
  ylab("Normalised log transformed counts")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size=18)) +
  scale_x_discrete(breaks=c("BS Nfree","HE Nfree","SS Nfree",
                            "BS Nitrogen","HE Nitrogen","SS Nitrogen"),
                   labels=c("BS Nfree","HE Nfree","SS Nfree",
                            "BS Nitrogen","HE Nitrogen","SS Nitrogen")) +
  facet_grid(Source~OTU,scales="free")

saveotu.sig$Factor3<- factor(saveotu.sig$Factor3,levels=c("BS Nfree","HE Nfree","SS Nfree",
                                                          "BS Nitrogen","HE Nitrogen","SS Nitrogen"))

pdf("xiaonaDEphase1.pdf")
for(i in significantotus){
  otusel <- saveotu.sig[saveotu.sig$OTU==i,]
  plotsig.vals <- ggplot(aes(y=count,x=Factor3,color=Factor3,shape=Factor3),
                         data=otusel) + 
    geom_boxplot(alpha=0.5,outlier.shape = 8,outlier.size = 3,
                 outlier.alpha = 1) +
    geom_jitter(position=position_jitter(width=.3,
                                         height=0),
                size=4,alpha=1) +
    scale_y_continuous(labels=scientific,trans='log10',
                       breaks = trans_breaks('log10',
                                             function(x) 10^x)) +
    ylab("Normalised log transformed counts")+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          text = element_text(size=18)) +
    scale_x_discrete(breaks=c("BS Nfree","HE Nfree","SS Nfree",
                              "BS Nitrogen","HE Nitrogen","SS Nitrogen"),
                     labels=c("BS Nfree","HE Nfree","SS Nfree",
                              "BS Nitrogen","HE Nitrogen","SS Nitrogen")) +
    facet_grid(~OTU,scales="free")  
  print(plotsig.vals)
  
}
dev.off()

#### 2. Phase 2 ####
meta <- metadata.endpoint
shared <- shared.endpoint
sh.x <- shared
sh.x.ns <- sh.x[which(rowSums(sh.x)!=1),]
meta <- data.frame(meta,Factor3=factor(paste(meta$Factor1,meta$Factor2))) # to have a one-factor DESeq design
## prevalence = (fraction of samples in which an OTU is observed minimum 1 time)
minobs=1
prevalence <- apply(as.matrix(sh.x.ns),1,function(x,minobs){sum(x>=minobs)},minobs)/ncol(sh.x.ns)
prevalencefilter <- prevalence>0.05
sh.x.prevfilt <- sh.x.ns[prevalencefilter,]
##Read counts should exceed 0.5 times the number of samples
sh.x.prevfilt <- sh.x.prevfilt[rowSums(sh.x.prevfilt)>0.5*ncol(sh.x.prevfilt),]
sh.x.prevfilt.prop <- sweep(sh.x.prevfilt,2,colSums(sh.x.prevfilt),'/')
deseqdata <- as.matrix(sh.x.prevfilt +1)
deseqdata <- DESeqDataSetFromMatrix(deseqdata,colData=meta,design= ~ Factor3)
deseqdata <- estimateSizeFactors(deseqdata)
sharedfiltereddeseq <- counts(deseqdata,normalized=TRUE)
deseqdata <- estimateDispersions(deseqdata,fitType="local")
sharedfiltereddeseqvartransf <- varianceStabilizingTransformation(deseqdata, blind = FALSE)
sharedfiltereddeseqvartransfmatrix <- assay(sharedfiltereddeseqvartransf)
wnwndeseqvartransf <- getVarianceStabilizedData(deseqdata)
testdiff <- DESeq(deseqdata,test = "Wald",fitType = "local")
testdiff.res5 <- results(testdiff,alpha=0.05,pAdjustMethod = "BY")
testdiff.res.ordered5 <- testdiff.res5[order(testdiff.res5$padj),]
summary(testdiff.res5)
signotu5 <- testdiff.res5$padj<0.05
resdf5 <- data.frame(sh.x.prevfilt,sigdif=signotu5,padjusted=testdiff.res5$padj)
write.csv(resdf5,"differentialxiaona5_endpoint.csv")
pwvals <- results(testdiff,contrast = c("Factor3","BS Nfree","BS Nitrogen"),pAdjustMethod = "BY")$padj
saveotu.sig <- data.frame()
significantotus <- rownames(resdf5[!is.na(resdf5$sigdif)&resdf5$sigdif==TRUE,])
for(i in significantotus){
  diffotu <- plotCounts(testdiff,gene=i,intgroup = "Factor3",returnData = TRUE,transform=FALSE)
  diffotu <- cbind(diffotu,OTU=i)
  saveotu.sig <- rbind(saveotu.sig,diffotu)
}
saveotu.sig <- data.frame(saveotu.sig,Source=factor(substr(as.character(saveotu.sig$Factor3),0,2)))
rankedsig <- saveotu.sig %>%  group_by(OTU) %>% summarise(med=median(count)) %>% arrange(desc(med)) %>% select(OTU)
plotsig.vals <- ggplot(aes(y=count,x=Factor3,color=Factor3,shape=Factor3),
                       data=saveotu.sig) + 
  geom_boxplot(alpha=0.5,outlier.shape = 8,outlier.size = 3,
               outlier.alpha = 1) +
  geom_jitter(position=position_jitter(width=.3,
                                       height=0),
              size=4,alpha=1) +
  scale_y_continuous(labels=scientific,trans='log10',
                     breaks = trans_breaks('log10',
                                           function(x) 10^x)) +
  ylab("Normalised log transformed counts")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size=18)) +
  scale_x_discrete(breaks=c("BS Nfree","HE Nfree","SS Nfree",
                            "BS Nitrogen","HE Nitrogen","SS Nitrogen"),
                   labels=c("BS Nfree","HE Nfree","SS Nfree",
                            "BS Nitrogen","HE Nitrogen","SS Nitrogen")) +
  facet_grid(Source~OTU,scales="free")

saveotu.sig$Factor3<- factor(saveotu.sig$Factor3,levels=c("BS Nfree","HE Nfree","SS Nfree",
                                                          "BS Nitrogen","HE Nitrogen","SS Nitrogen"))

pdf("xiaonaDEendpoint.pdf")
for(i in significantotus){
  otusel <- saveotu.sig[saveotu.sig$OTU==i,]
  plotsig.vals <- ggplot(aes(y=count,x=Factor3,color=Factor3,shape=Factor3),
                         data=otusel) + 
    geom_boxplot(alpha=0.5,outlier.shape = 8,outlier.size = 3,
                 outlier.alpha = 1) +
    geom_jitter(position=position_jitter(width=.3,
                                         height=0),
                size=4,alpha=1) +
    scale_y_continuous(labels=scientific,trans='log10',
                       breaks = trans_breaks('log10',
                                             function(x) 10^x)) +
    ylab("Normalised log transformed counts")+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          text = element_text(size=18)) +
    scale_x_discrete(breaks=c("BS Nfree","HE Nfree","SS Nfree",
                              "BS Nitrogen","HE Nitrogen","SS Nitrogen"),
                     labels=c("BS Nfree","HE Nfree","SS Nfree",
                              "BS Nitrogen","HE Nitrogen","SS Nitrogen")) +
    facet_grid(~OTU,scales="free")  
  print(plotsig.vals)
  
}
dev.off()