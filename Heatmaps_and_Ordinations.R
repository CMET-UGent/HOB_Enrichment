####                   Heatmaps_and_Ordinations.R                          #####
#      version: 2019-11-01, created by: Frederiek - Maarten Kerckhof           #
#      Analysis code for Hu et al. (2019), released under GPL-3                #
################################################################################


#### load required packages ---------------------------------------------------- 
library(ggplot2)
library(extrafont)
library(phyloseq)
library(NMF)
library(vegan)
library(gplots)   
library(readxl)
library(DESeq2)
library(data.table)
library(splitstackshape)
library(dplyr)
library(scales)
library(ape)
#### custom functions ----------------------------------------------------------
preformattax <- function(taxonomy)
{
  #Input: mothur taxonomy with probs=T (default)
  #Output: splitted and stripped tax file for R processing (e.g. phyloseq)
  tax.good <- cSplit(taxonomy,"Taxonomy",";")
  tax.good.probs <- do.call(cbind, lapply(tax.good[,3:ncol(tax.good)], 
                                          function(x){
                                            do.call(rbind, 
                                                    strsplit(as.character(x), 
                                                             "\\((?=\\d)", 
                                                             perl = TRUE))
                                          }
  )
  )
  # tax.good.probs.nona <- subset(tax.good.probs,
  #           select=which(colSums(apply(tax.good.probs,2,is.na))==0))
  tax.final <- as.data.frame(apply(tax.good.probs,2,
                                   function(x) sub(")","",x)))
  otunames <- tax.good$OTU
  rownames(tax.final) <- otunames
  colnames(tax.final) <- c("Regnum","Prob_R","Phylum","Prob_P",
                           "Classis","Prob_C","Ordo","Prob_O",
                           "Familia","Prob_F","Genus","Prob_G")
  if(ncol(tax.final)==12){
    #for RDP & SILVA taxonomies which stop at genus level
    tax.final <- as.data.frame(mutate(tax.final,Species=NA,Prob_S=NA))
    rownames(tax.final) <- otunames
  }else{
    #for greengenes taxonomy
    colnames(tax.final)[13:14] <- c("Species","Prob_S")
    tax.final <- as.data.frame(apply(tax.final,
                                     2,
                                     function(x)sub("[a-z]__","",x)))
    tax.final$Species <- paste(tax.final$Genus,tax.final$Species) 
    #just the species is not very sensible
  }
  #to comply to mothur 1.38 tax format
  charvectgenus <- as.character(tax.final$Genus)
  probvectgenus <- as.numeric(as.character(tax.final$Prob_G))
  nasgenuslev <- which(is.na(charvectgenus))
  for(i in nasgenuslev)
  {
    if(grepl("unclassified",as.character(tax.final$Familia[i])))
    {
      if(!is.na(as.character(tax.final$Familia[i])))
      {
        charvectgenus[i]<-as.character(tax.final$Familia[i])
      }else{
        if(grepl("unclassified",as.character(tax.final$Ordo[i])))
        {
          charvectgenus[i]<-as.character(tax.final$Ordo[i])
        }
      }
      
    }else{
      if(!is.na(as.character(tax.final$Familia[i])))
      {
        charvectgenus[i]<-paste0(as.character(tax.final$Familia[i]),"_unclassified")
      }else{
        charvectgenus[i]<-paste0(as.character(tax.final$Ordo[i]),"_unclassified")
      }
      
    }
    probvectgenus[i] <- as.numeric(as.character(tax.final$Prob_F[i]))
  }
  tax.final$Genus <- factor(charvectgenus)
  tax.final$Prob_G <- probvectgenus
  
  charvectfamilia <- as.character(tax.final$Familia)
  probvectfamilia <- as.numeric(as.character(tax.final$Prob_F))
  nasfamlev <- which(is.na(charvectfamilia))
  for(i in nasfamlev)
  {
    if(grepl("unclassified",as.character(tax.final$Ordo[i])))
    {
      if(!is.na(as.character(tax.final$Ordo[i])))
      {
        charvectfamilia[i]<-as.character(tax.final$Ordo[i])
      }else{
        if(grepl("unclassified",as.character(tax.final$Classis[i])))
        {
          charvectfamilia[i]<-as.character(tax.final$Classis[i])
        }
      }
      
    }else{
      if(!is.na(as.character(tax.final$Ordo[i])))
      {
        charvectfamilia[i]<-paste0(as.character(tax.final$Ordo[i]),"_unclassified")
      }else{
        charvectfamilia[i]<-paste0(as.character(tax.final$Classis[i]),"_unclassified")
      }
      
    }
    probvectfamilia[i] <- as.numeric(as.character(tax.final$Prob_O[i]))
  }
  tax.final$Familia <- factor(charvectfamilia)
  tax.final$Prob_F <- probvectfamilia
  return(tax.final)
}


#### load the data -------------------------------------------------------------
# Read shared file (OTU table)
fldr <- "Source_Data"
fn <- "stability"
shared <- fread(input = "Source_Data/HOB_final_shared.shared",
                header=TRUE) 
#replace here with the filename and location of your mothur shared file
shared <- as.data.frame(shared) 
# fread reads in as tibble, convert to data.frame for downstream processing
desgroups <- shared$Group 
# select sample names
shared.x <- shared[,4:ncol(shared)] 
# remove all non-count data from the OTU table
rownames(shared.x) <- desgroups
shared.t <- as.data.frame(t(shared.x))
# Removing absolute singletons
shared.t.ns <- shared.t[which(rowSums(shared.t)!=1),]
sharedns <- data.frame(label=rep(0.03,ncol(shared.t.ns)),
                       Group=colnames(shared.t.ns),
                       numOtus=nrow(shared.t.ns),t(shared.t.ns))
# read taxonomy file
otutaxonomy <- fread(input = "Source_Data/HOB_final_taxonomy.taxonomy",
                     header=TRUE) 
# replace here with the filename and location of your mothur taxonomy file

taxonomy.spl <- preformattax(otutaxonomy)
# Remove absolute singleton OTUs from the taxonomy as well
taxonomy.np <- taxonomy.spl %>% dplyr::select(-dplyr::contains("Prob"))

#filter taxonomy
taxonomy.np.ns <- taxonomy.np[which(rownames(taxonomy.np) 
                                    %in% rownames(shared.t.ns)),]
metadata.tibble <- readxl::read_excel("Source_Data/MetaData_5.xlsx",
                                      sheet="ForR") 
#set location of your own metadata.xlsx, which should contain a sheet "ForR"
factdescs <- readxl::read_excel("Source_Data/MetaData_5.xlsx",
                                sheet="FactDesc")
metadata <- as.data.frame(metadata.tibble) 
#to avoid warnings/errors with rownames
if(is.numeric(metadata$SampleName)) #check for fully numeric sample names
{
  metadata$SampleName <- paste(dataname,metadata$SampleName,sep="")
}
rownames(metadata) <- metadata$SampleName

colnames(shared.t.ns) <- plyr::mapvalues(colnames(shared.t.ns),
                                         from=as.character(metadata$Code),
                                         to=as.character(metadata$SampleName)) 
#rename

shared.t.ns <- shared.t.ns %>% select(-c(`SS0-a`,`SS0-b`,`SS0-c`))
metadata <- metadata %>%  filter(!SampleName%in%c("SS0-a","SS0-b","SS0-c"))
rownames(metadata) <- metadata$SampleName
metadata.smpdat <- sample_data(metadata)

# read celcount info
cellconcentations <- read.csv("Source_Data/Cellcounts.csv",header = FALSE, 
                              col.names = c("Sample","ConcentrationpermL"))

### phyloseq constructors ------------------------------------------------------
otumat.ns <- as.matrix(shared.t.ns)
taxmat.ns <- as.matrix(taxonomy.np.ns)
OTU       <- otu_table(otumat.ns,taxa_are_rows = TRUE)
TAX       <- tax_table(taxmat.ns)
physeqobj <- phyloseq(OTU,TAX)
physeqobj.meta <- merge_phyloseq(physeqobj,metadata.smpdat)

### prevalence filtering according to supplementary of McMurdie & Holmes (2104)-
## prevalence = (fraction of samples in which an OTU is observed minimum 1 time)
minobs=1
prevalence <- apply(as.matrix(shared.t.ns),1,
                    function(x,minobs){sum(x>=minobs)},minobs)/ncol(shared.t.ns)
prevalencefilter <- prevalence>0.05
sharedminsingletonwnwn <- shared.t.ns[prevalencefilter,]
##Read counts should exceed 0.5 times the number of samples
sharedfilteredwnwn <- 
  sharedminsingletonwnwn[rowSums(sharedminsingletonwnwn)>0.5*
                           ncol(sharedminsingletonwnwn),]
# deseq normalise ==> very dependent upon design
metadata$factor1 <- factor(metadata$Factor1)
metadata$factor2 <- factor(metadata$Factor2)
deseqdata <- as.matrix(sharedfilteredwnwn +1)

#### heatmaps without standardization ------------------------------------------
relabsphyseq <- transform_sample_counts(physeqobj, function(x) x/sum(x)) 
# create relative abudances
Top24OTUs.relab <- names(sort(taxa_sums(relabsphyseq), TRUE)[1:24])      
# select top 24 OTUs in relative abundance
physeqobj24.relab  <- prune_taxa(Top24OTUs.relab, relabsphyseq)          
# extract the top 24 OTUs
relab24mat <- as.matrix(otu_table(physeqobj24.relab))                    
# Create a matrix with the selected relative abundances
rownames(relab24mat) <- 
  as.character(as.data.frame(tax_table(physeqobj24.relab)@.Data)$Genus)
my_palette <- colorRampPalette(c("black","green" ,"red"))(n = 409)
# this colors can be selected from a hexadecimal color picker
# e.g. https://www.google.be/search?hl=en&q=hex+color+picker&gws_rd=ssl

#heatmap.2(relab24mat,col=brewer.pal(9,name="Reds"))

col_breaks = c(seq(0,0.050,length=200),seq(0.051,0.10,length=100),
               seq(0.11,0.15,length=50),seq(0.151,0.20,length=50),
               seq(0.21,max(relab24mat),length=10))

heatmap.2(relab24mat, main = "Genus level heatmap",
          notecol="black", density.info="none",
          trace="none", margins =c(10,10), col=my_palette, 
          dendrogram="column", breaks=col_breaks,
          lmat=rbind(c(3,4),c(1,2)),lwi = c(3,1),
          Rowv = FALSE, Colv = TRUE)

#### family-level heatmap ------------------------------------------------------
relab24matfam <- as.matrix(otu_table(physeqobj24.relab))  
rownames(relab24matfam) <- 
 paste(as.character(as.data.frame(unclass(tax_table(physeqobj24.relab)))$Familia),
       as.character(as.data.frame(unclass(tax_table(physeqobj24.relab)))$Genus))


col_breaksfam = c(seq(0,0.050,length=200),seq(0.051,0.10,length=100),
                  seq(0.11,0.15,length=50),seq(0.151,0.20,length=50),
                  seq(0.21,max(relab24matfam),length=10))

heatmap.2(relab24matfam, main = "Family and Genus level heatmap",
          notecol="black", density.info="none",
          trace="none", margins =c(10,10), col=my_palette, 
          dendrogram="column", breaks=col_breaksfam,
          lmat=rbind(c(3,4),c(1,2)),lwi = c(3,1),
          Rowv = FALSE, Colv = TRUE)

#### top 100 more advanced heatmap without standardization ---------------------
physeqobj.meta.subs <- prune_taxa(names(sort(taxa_sums(physeqobj.meta),
                                             TRUE)[1:10]),
                                  physeqobj.meta)
plot_heatmap(physeqobj.meta.subs, sample.label="SampleName",
             taxa.label="Genus",method = "NMDS",distance = "jaccard",
             low="#294ddb", high="#db2929", na.value="white",
             trans = scales::log_trans(10))


#### heatmaps with standardization ---------------------------------------------
topn <- 20
wnwnsh <- sharedfilteredwnwn
wnwnsh.prop <- decostand(wnwnsh,method = "total",MARGIN=2)
shns <- shared.t.ns
shns.prop <- decostand(shns,method = "total",MARGIN=2)
shns.log <- decostand(shns,method = "log",logbase=10)
shns.cs <- shns.prop*as.numeric(quantile(colSums(shared.t.ns),.75)) 
# other quantile then max see: #boxplot(colSums(shared.t.ns))
wnwnsh.log <- decostand(wnwnsh,method = "log",logbase=10)

metadata.colsidecols <- metadata %>% select(c(Factor1,Factor2)) %>%
  transmute(Source=Factor1,
            Nitrogen=as.factor(Factor2))
rownames(metadata.colsidecols) <- metadata$SampleName
md.colsidecols.ordered <- metadata.colsidecols[colnames(wnwnsh),]

wnwnsh.top <- wnwnsh %>% dplyr::mutate(otu=rownames(wnwnsh)) %>% 
  top_n(topn,rowSums(wnwnsh[1:(ncol(wnwnsh)-1)]))
rownames(wnwnsh.top) <- wnwnsh.top$otu
wnwnsh.top <- wnwnsh.top %>% dplyr::select(-c(otu))
wnwnsh.top.clust <- hclust(vegdist(t(wnwnsh.top),"jaccard",binary = FALSE),
                           method="ward.D2")

NMF::aheatmap(wnwnsh.top,border_color="white",Rowv=NA,
              Colv=wnwnsh.top.clust$order,
              annCol=md.colsidecols.ordered, 
              main="WNWN-filtered abundance-based jaccard Ward D2")

wnwnsh.prop.top <- wnwnsh.prop %>% dplyr::mutate(otu=rownames(wnwnsh.prop)) %>% 
  top_n(topn,rowSums(wnwnsh.prop[1:(ncol(wnwnsh.prop)-1)]))
rownames(wnwnsh.prop.top) <- wnwnsh.prop.top$otu
wnwnsh.prop.top <- wnwnsh.prop.top %>% dplyr::select(-c(otu))
wnwnsh.prop.top.clust <- hclust(vegdist(t(wnwnsh.prop.top),"jaccard",
                                        binary = FALSE),method="ward.D2")

NMF::aheatmap(wnwnsh.prop.top,border_color="white",Rowv=NA,
              Colv=wnwnsh.prop.top.clust$order,
              annCol=md.colsidecols.ordered,
              main="WNWN-filtered (abun) jaccard Ward D2 relative abundances")

#### NMDS and PCoA -------------------------------------------------------------
set.seed(4798) # for reproducible randomness
# first, rescale the data as appropriate
shared.relab <- decostand(shared.t.ns,method="total",MARGIN=2)
shared.cs <- shared.relab*min(colSums(shared.t.ns)) #common scale strategy
# then, we calculate the distances on the rescaled counts
disjac <- vegdist(t(shared.relab),method="jaccard") 
#calculate distances between samples
disjac.cs <- vegdist(t(shared.t.ns),method="jaccard")
discyd <- vegdist(t(shared.t.ns),method="cao") 
#can handle different sample sizes
suppressWarnings(discyd.cs <- vegdist(t(ceiling(shared.cs)),method="cao")) 
#can handle different sample sizes

# Calculate absolute counts based upon flow cytometry
shared.absolute <- round(t(t(shared.relab)*
                             cellconcentations$ConcentrationpermL))

# Run the ordination methods (NMDS ) on the data without absolute counts
mdsjac <- metaMDS(disjac,trace=FALSE,try=50,trymax=2000)
mdsjac.cs <- metaMDS(disjac.cs,trace=FALSE,try=50,trymax=5000)
mdscyd <- metaMDS(discyd,trace=FALSE,try=50,trymax=5000)
mdscyd.cs <- metaMDS(discyd.cs,trace=FALSE,try=50,trymax=2000)

# then, we calculate the distances on the absolute counts
disjac.abs <- vegdist(t(shared.absolute),method="jaccard")
discyd.abs <- vegdist(t(shared.absolute),method="cao") 
#can handle different sample sizes
# Run the ordination methods (NMDS ) on the data WITH absolute counts
mdsjac.abs <- metaMDS(disjac.abs,trace=FALSE,try=50,trymax=2000)
mdscyd.abs <- metaMDS(discyd.abs,trace=FALSE,try=50,trymax=5000)


# here, we plot the output for the nmds on the absolute count data
plot(mdsjac.abs,type="n",display="sites",
     main="NMDS, Jaccard distance, absolute counts")
ordipointlabel(mdsjac.abs,display = "sites",add=TRUE)
colvect<-rep("blue",length(rownames(mdsjac$points)))
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[2])]<-"red"
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[3])]<-"green"
points(mdsjac.abs,pch=16,col=colvect)
legend("bottomright",legend = levels(factor(metadata$Factor2)),
       col=c("blue","red","green"),pch = 16)



plot(mdscyd.abs,type="n",display="sites",
     main="NMDS, Cao distance, absolute counts")
ordipointlabel(mdscyd.abs,display = "sites",add=TRUE)
colvect<-rep("blue",length(rownames(mdscyd.abs$points)))
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[2])]<-"red"
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[3])]<-"green"
points(mdscyd.abs,pch=16,col=colvect)
legend("bottomright",legend = levels(factor(metadata$Factor2)),
       col=c("blue","red","green"),pch = 16)


# run the pcoa method
# based upon the jaccard distance and the relative abundance data, 
# with jaccard distance
pcoajac <- pcoa(disjac,correction="cailliez")
biplot(pcoajac,t(shared.relab))
# based upon the jaccard distance and the absolute abundance data
pcoajac.abs <- pcoa(disjac.abs,correction="cailliez")
biplot(pcoajac,t(shared.absolute))

# Pcoa through vegan cmdscale (Euclidean)
# based upon the relative abundance data, with jaccard distance
pcoajac <- cmdscale(disjac, eig=TRUE,x.ret=TRUE,list.=TRUE)
ordiplot(pcoajac,type="n",display="sites", 
         main="PCoa based on relative abundances and Jaccard distance")
ordipointlabel(pcoajac,display = "sites",add=TRUE)
colvect<-rep("blue",length(rownames(pcoajac$points)))
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[2])]<-"red"
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[3])]<-"green"
points(pcoajac$points,pch=16,col=colvect)
legend("bottomright",legend = levels(factor(metadata$Factor2)),
       col=c("blue","red","green"),pch = 16)

pcoacyd <- cmdscale(discyd, eig=TRUE,x.ret=TRUE,list.=TRUE)
ordiplot(pcoacyd,type="n",display="sites", 
         main="PCoa based on relative abundances and Cao distance")
ordipointlabel(pcoacyd,display = "sites",add=TRUE)
colvect<-rep("blue",length(rownames(pcoacyd$points)))
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[2])]<-"red"
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[3])]<-"green"
points(pcoacyd$points,pch=16,col=colvect)
legend("bottomright",legend = levels(factor(metadata$Factor2)),
       col=c("blue","red","green"),pch = 16)

# based upon absolute abundance data, with jaccard distance
pcoajac.abs <- cmdscale(disjac.abs)
plot(pcoajac.abs,type="n")
ordipointlabel(pcoajac.abs,display = "sites",add=TRUE)
colvect<-rep("blue",length(rownames(mdsjac.abs$points)))
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[2])]<-"red"
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[3])]<-"green"
points(pcoajac.abs,pch=16,col=colvect)
legend("topright",legend = levels(factor(metadata$Factor2)),
       col=c("blue","red","green"),pch = 16)

ordiplot(pcoajac.abs,type = "text",display="sites")

biplot(pcoajac.abs,t(shared.absolute))

#### NMDS and PCoA without inoculum --------------------------------------------
set.seed(4798)

shared.no.inoculum <- shared.t.ns %>% dplyr::select(-contains("0"))
shared.no.inoculum.abs <- as.data.frame(shared.absolute) %>% 
  dplyr::select(-contains("0"))
metadata.no.inoculum <- metadata %>% 
  dplyr::filter(!grepl("0",rownames(metadata)))
rownames(metadata.no.inoculum) <- 
  rownames(metadata)[!grepl("0",rownames(metadata))]



shared.relab <- decostand(shared.no.inoculum,method="total",MARGIN=2) 
# make relative abundances
shared.cs <- shared.relab*min(colSums(shared.t.ns)) # make common scale
disjac <- vegdist(t(shared.relab),method="jaccard")
disjac.cs <- vegdist(t(shared.t.ns),method="jaccard")
discyd <- vegdist(t(shared.t.ns),method="cao") 
#can handle different sample sizes
suppressWarnings(discyd.cs <- vegdist(t(ceiling(shared.cs)),method="cao")) 
#can handle different sample sizes

mdsjac <- metaMDS(disjac,trace=FALSE,try=50,trymax=2000)
mdsjac.cs <- metaMDS(disjac.cs,trace=FALSE,try=50,trymax=5000)
mdscyd <- metaMDS(discyd,trace=FALSE,try=50,trymax=5000)
mdscyd.cs <- metaMDS(discyd.cs,trace=FALSE,try=50,trymax=2000)


plot(mdsjac,type="n",display="sites")
ordipointlabel(mdsjac,display = "sites",add=TRUE)

colvect <- rep("blue",length(rownames(mdsjac$points)))
colvect[which(metadata.no.inoculum$Factor2==
                levels(factor(metadata.no.inoculum$Factor2))[2])]<-"red"


pchvect <- rep(15,length(rownames(mdsjac$points)))
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[2])]<-16
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[3])]<-17

points(mdsjac,pch=pchvect,col=colvect,cex=1.5)
legend("top",
       legend = paste(rep(levels(factor(metadata.no.inoculum$Factor1)),2),
                      levels(factor(metadata.no.inoculum$Factor2))),
       col=rep(c("blue","red"),3),pch = rep(15:17,2))

plot(mdsjac.cs,type="n",display="sites")
ordipointlabel(mdsjac.cs,display = "sites",add=TRUE)
colvect<-rep("blue",length(rownames(mdsjac.cs$points)))
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[2])]<-"red"
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[3])]<-"green"
points(mdsjac.cs,pch=16,col=colvect)
legend("bottomright",legend = levels(factor(metadata$Factor2)),
       col=c("blue","red","green"),pch = 16)

plot(mdscyd,type="n",display="sites")
ordipointlabel(mdscyd,display = "sites",add=TRUE)
colvect<-rep("blue",length(rownames(mdscyd$points)))
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[2])]<-"red"
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[3])]<-"green"
points(mdscyd,pch=16,col=colvect)
legend("bottomright",legend = levels(factor(metadata$Factor2)),
       col=c("blue","red","green"),pch = 16)

plot(mdscyd.cs,type="n",display="sites")
ordipointlabel(mdscyd.cs,display = "sites",add=TRUE)
colvect<-rep("blue",length(rownames(mdscyd.cs$points)))
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[2])]<-"red"
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[3])]<-"green"
points(mdscyd.cs,pch=16,col=colvect)
legend("bottomright",legend = levels(factor(metadata$Factor2)),
       col=c("blue","red","green"),pch = 16)


stressplot(mdsjac,disjac)
stressplot(mdsjac.cs,disjac.cs)
stressplot(mdscyd,discyd)
stressplot(mdscyd.cs,discyd.cs)


pcoajac <- pcoa(disjac,correction="cailliez")
biplot(pcoajac,t(shared.relab))

pcoajac <- cmdscale(disjac)
plot(pcoajac,type="n")
ordipointlabel(pcoajac,display = "sites",add=TRUE)
colvect<-rep("blue",length(rownames(pcoajac)))
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[2])]<-"red"
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[3])]<-"green"
points(pcoajac,pch=16,col=colvect)
legend("bottomright",legend = levels(factor(metadata$Factor2)),col=c("blue","red","green"),pch = 16)

#### no inoculum absolute counts -----------------------------------------------
disjac.abs.ni <- vegdist(t(shared.no.inoculum.abs),method="jaccard")
discyd.abs.ni <- vegdist(t(shared.no.inoculum.abs),method="cao") 
#can handle different sample sizes
mdsjac.abs.ni <- metaMDS(disjac.abs.ni,trace=FALSE,try=50,trymax=2000)
mdscyd.abs.ni <- metaMDS(discyd.abs.ni,trace=FALSE,try=50,trymax=5000)

plot(mdsjac.abs.ni,type="n",display="sites")
ordipointlabel(mdsjac.abs.ni,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdsjac.abs.ni$points)))
colvect[which(metadata.no.inoculum$Factor2==
                levels(factor(metadata.no.inoculum$Factor2))[2])]<-"red"
pchvect <- rep(15,length(rownames(mdsjac.abs.ni$points)))
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[2])]<-16
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[3])]<-17
points(mdsjac.abs.ni,pch=pchvect,col=colvect,cex=1.5)
legend("top",
       legend = paste(rep(levels(factor(metadata.no.inoculum$Factor1)),2),
                      levels(factor(metadata.no.inoculum$Factor2))),
       col=rep(c("blue","red"),3),pch = rep(15:17,2))

plot(mdscyd.abs.ni,type="n",display="sites")
ordipointlabel(mdscyd.abs.ni,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdscyd.abs.ni$points)))
colvect[which(metadata.no.inoculum$Factor2==
                levels(factor(metadata.no.inoculum$Factor2))[2])]<-"red"
pchvect <- rep(15,length(rownames(mdscyd.abs.ni$points)))
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[2])]<-16
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[3])]<-17
points(mdscyd.abs.ni,pch=16,col=colvect)
legend("topright",
       legend = paste(rep(levels(factor(metadata.no.inoculum$Factor1)),2),
                      levels(factor(metadata.no.inoculum$Factor2))),
       col=rep(c("blue","red"),3),pch = rep(15:17,2))


pcoajac <- pcoa(disjac.abs.ni,correction="cailliez")
biplot(pcoajac,t(shared.no.inoculum.abs))


pcoajac.abs.ni <- cmdscale(disjac.abs.ni,eig=TRUE)
eigvalpcoajac.abs.ni <- scale(abs(pcoajac.abs.ni$eig),
                              center=FALSE,scale=sum(abs(pcoajac.abs.ni$eig)))
plot(pcoajac.abs.ni$points,type="n",
     xlab=paste("PCoA axis 1 (",signif(eigvalpcoajac.abs.ni[1]*100,4),"%)"),
     ylab=paste("PCoA axis 2 (",signif(eigvalpcoajac.abs.ni[2]*100,4),"%)"))
#æordipointlabel(pcoajac.abs.ni,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdscyd.abs.ni$points)))
colvect[which(metadata.no.inoculum$Factor2==
                levels(factor(metadata.no.inoculum$Factor2))[2])]<-"red"
pchvect <- rep(15,length(rownames(mdscyd.abs.ni$points)))
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[2])]<-16
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[3])]<-17
points(pcoajac.abs.ni$points,pch=pchvect,col=alpha(colvect,0.6))

legend("topleft",
       legend= expression(paste("HE-",N[2]),paste("BS-",N[2]),
                          paste("SS-",N[2]),paste("HE-",NH[4]^"+"),
                          paste("BS-",NH[4]^"+"),paste("SS-",NH[4]^"+")),
       col=alpha(c("blue","blue","blue","red","red","red"),0.6),
       pch = c(16,15,17,16,15,17),
       bty="n")


pcoajac.abs.ni <- cmdscale(disjac.abs.ni,eig=TRUE)
eigvalpcoajac.abs.ni <- scale(abs(pcoajac.abs.ni$eig),
                              center=FALSE,scale=sum(abs(pcoajac.abs.ni$eig)))
plot(pcoajac.abs.ni$points,type="n",
     xlab=paste("PCoA axis 1 (",signif(eigvalpcoajac.abs.ni[1]*100,4),"%)"),
     ylab=paste("PCoA axis 2 (",signif(eigvalpcoajac.abs.ni[2]*100,4),"%)"))
#ordipointlabel(pcoajac.abs.ni,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdscyd.abs.ni$points)))
colvect[which(metadata.no.inoculum$Factor2==
                levels(factor(metadata.no.inoculum$Factor2))[2])]<-"red"
pchvect <- rep(15,length(rownames(mdscyd.abs.ni$points)))
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[1]&
                grepl("1",metadata.no.inoculum$SampleName))]<-0
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[2]&
                grepl("2",metadata.no.inoculum$SampleName))]<-16
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[2]&
                grepl("1",metadata.no.inoculum$SampleName))]<-1
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[3]&
                grepl("2",metadata.no.inoculum$SampleName))]<-17
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[3]&
                grepl("1",metadata.no.inoculum$SampleName))]<-2
points(pcoajac.abs.ni$points,pch=pchvect,col=alpha(colvect,0.6),cex=2)

# legend("topleft",
#        legend= expression(paste("HE-",N[2]),paste("BS-",N[2]),
#                           paste("SS-",N[2]),paste("HE-",NH[4]^"+"),
#                           paste("BS-",NH[4]^"+"),paste("SS-",NH[4]^"+")),
#        col=alpha(c("blue","blue","blue","red","red","red"),0.6),
#        pch = c(16,15,17,16,15,17),
#        bty="n")





tiff(filename = "PCoA_jaccard_abs_noInoculum_3in.tiff",width = 3,height = 3,
     units = "in",pointsize = 8,family = "sans",res=600)
plot(pcoajac.abs.ni$points,type="n",
     xlab=paste("PCoA axis 1 (",signif(eigvalpcoajac.abs.ni[1]*100,4),"%)"),
     ylab=paste("PCoA axis 2 (",signif(eigvalpcoajac.abs.ni[2]*100,4),"%)"))
#ordipointlabel(pcoajac.abs.ni,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdscyd.abs.ni$points)))
colvect[which(metadata.no.inoculum$Factor2==
                levels(factor(metadata.no.inoculum$Factor2))[2])]<-"red"
pchvect <- rep(15,length(rownames(mdscyd.abs.ni$points)))
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[1]&
                grepl("1",metadata.no.inoculum$SampleName))]<-0
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[2]&
                grepl("2",metadata.no.inoculum$SampleName))]<-16
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[2]&
                grepl("1",metadata.no.inoculum$SampleName))]<-1
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[3]&
                grepl("2",metadata.no.inoculum$SampleName))]<-17
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[3]&
                grepl("1",metadata.no.inoculum$SampleName))]<-2
points(pcoajac.abs.ni$points,pch=pchvect,col=alpha(colvect,0.6),cex=1)
dev.off()

tiff(filename = "PCoA_jaccard_abs_noInoculum_LZW_compressed_3in.tiff",width = 3,
     height = 3,units = "in",pointsize = 8,family = "sans",res=600,
     compression = "lzw" )
plot(pcoajac.abs.ni$points,type="n",
     xlab=paste("PCoA axis 1 (",signif(eigvalpcoajac.abs.ni[1]*100,4),"%)"),
     ylab=paste("PCoA axis 2 (",signif(eigvalpcoajac.abs.ni[2]*100,4),"%)"))
#ordipointlabel(pcoajac.abs.ni,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdscyd.abs.ni$points)))
colvect[which(metadata.no.inoculum$Factor2==
                levels(factor(metadata.no.inoculum$Factor2))[2])]<-"red"
pchvect <- rep(15,length(rownames(mdscyd.abs.ni$points)))
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[1]&
                grepl("1",metadata.no.inoculum$SampleName))]<-0
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[2]&
                grepl("2",metadata.no.inoculum$SampleName))]<-16
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[2]&
                grepl("1",metadata.no.inoculum$SampleName))]<-1
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[3]&
                grepl("2",metadata.no.inoculum$SampleName))]<-17
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[3]&
                grepl("1",metadata.no.inoculum$SampleName))]<-2
points(pcoajac.abs.ni$points,pch=pchvect,col=alpha(colvect,0.6),cex=1)
dev.off()

tiff(filename = "PCoA_jaccard_abs_noInoculum_LZW_compressed_withlabels_3in.tiff",
     width = 3,height = 3,units = "in",pointsize = 8,family = "sans",res=600,
     compression = "lzw" )
plot(pcoajac.abs.ni$points,type="n",
     xlab=paste("PCoA axis 1 (",signif(eigvalpcoajac.abs.ni[1]*100,4),"%)"),
     ylab=paste("PCoA axis 2 (",signif(eigvalpcoajac.abs.ni[2]*100,4),"%)"))
ordipointlabel(pcoajac.abs.ni,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdscyd.abs.ni$points)))
colvect[which(metadata.no.inoculum$Factor2==
                levels(factor(metadata.no.inoculum$Factor2))[2])]<-"red"
pchvect <- rep(15,length(rownames(mdscyd.abs.ni$points)))
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[1]&
                grepl("1",metadata.no.inoculum$SampleName))]<-0
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[2]&
                grepl("2",metadata.no.inoculum$SampleName))]<-16
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[2]&
                grepl("1",metadata.no.inoculum$SampleName))]<-1
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[3]&
                grepl("2",metadata.no.inoculum$SampleName))]<-17
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[3]&
                grepl("1",metadata.no.inoculum$SampleName))]<-2
points(pcoajac.abs.ni$points,pch=pchvect,col=alpha(colvect,0.6),cex=1)
dev.off()


biplot(pcoajac)

pcoajac.abs.ni.wcm <- wcmdscale(disjac.abs.ni,eig=TRUE,k=2,add = "lingoes")
eigvalpcoajac.abs.ni.wcm <- scale(abs(pcoajac.abs.ni.wcm$eig),
                                  center=FALSE,
                                  scale=sum(abs(pcoajac.abs.ni.wcm$eig)))
a<-ordiplot(pcoajac.abs.ni.wcm,type="n",
            xlab=paste("PCoA Dim 1 (",
                       signif(eigvalpcoajac.abs.ni.wcm[1]*100,4),"%)",sep=""),
            ylab=paste("PCoA Dim 2 (",
                       signif(eigvalpcoajac.abs.ni.wcm[2]*100,4),"%)",sep=""))
ordipointlabel(pcoajac.abs.ni.wcm,"sites",add=TRUE)
colvect <- rep("blue",length(rownames(pcoajac.abs.ni.wcm$points)))
colvect[which(metadata.no.inoculum$Factor2==
                levels(factor(metadata.no.inoculum$Factor2))[2])]<-"red"
pchvect <- rep(15,length(rownames(pcoajac.abs.ni.wcm$points)))
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[2])]<-16
pchvect[which(metadata.no.inoculum$Factor1==
                levels(factor(metadata.no.inoculum$Factor1))[3])]<-17
points(pcoajac.abs.ni$points,pch=pchvect,col=colvect)
legend("top",
       legend = paste(rep(levels(factor(metadata.no.inoculum$Factor1)),2),
                      levels(factor(metadata.no.inoculum$Factor2))),
       col=rep(c("blue","red"),3),pch = rep(15:17,2))

#### endpoint NMDS and PCoA ----------------------------------------------------
set.seed(4798)

shared.endpoint <- shared.t.ns %>% 
  dplyr::select(which(!grepl(pattern = "[NES]1",colnames(shared.t.ns)))) %>%  
  dplyr::select(-contains("0"))
shared.endpoint.abs <- as.data.frame(shared.absolute) %>% 
  dplyr::select(which(!grepl(pattern = "[NES]1",colnames(shared.t.ns)))) %>%  
  dplyr::select(-contains("0"))
metadata.endpoint <- metadata[metadata$SampleName %in% colnames(shared.endpoint),]


shared.relab <- decostand(shared.endpoint,method="total",MARGIN=2) # make relative abundances
shared.cs <- shared.relab*min(colSums(shared.endpoint)) # make common scale
disjac <- vegdist(t(shared.relab),method="jaccard")
disjac.cs <- vegdist(t(shared.t.ns),method="jaccard")
discyd <- vegdist(t(shared.t.ns),method="cao") #can handle different sample sizes
suppressWarnings(discyd.cs <- vegdist(t(ceiling(shared.cs)),method="cao")) #can handle different sample sizes

mdsjac <- metaMDS(disjac,trace=FALSE,try=50,trymax=2000)
mdsjac.cs <- metaMDS(disjac.cs,trace=FALSE,try=50,trymax=5000)
mdscyd <- metaMDS(discyd,trace=FALSE,try=50,trymax=5000)
mdscyd.cs <- metaMDS(discyd.cs,trace=FALSE,try=50,trymax=2000)


plot(mdsjac,type="n",display="sites")
ordipointlabel(mdsjac,display = "sites",add=TRUE)

colvect <- rep("blue",length(rownames(mdsjac$points)))
colvect[which(metadata.no.inoculum$Factor2==levels(factor(metadata.no.inoculum$Factor2))[2])]<-"red"


pchvect <- rep(15,length(rownames(mdsjac$points)))
pchvect[which(metadata.no.inoculum$Factor1==levels(factor(metadata.no.inoculum$Factor1))[2])]<-16
pchvect[which(metadata.no.inoculum$Factor1==levels(factor(metadata.no.inoculum$Factor1))[3])]<-17

points(mdsjac,pch=pchvect,col=colvect,cex=1.5)
legend("top",
       legend = paste(rep(levels(factor(metadata.no.inoculum$Factor1)),2),
                      levels(factor(metadata.no.inoculum$Factor2))),
       col=rep(c("blue","red"),3),pch = rep(15:17,2))

plot(mdsjac.cs,type="n",display="sites")
ordipointlabel(mdsjac.cs,display = "sites",add=TRUE)
colvect<-rep("blue",length(rownames(mdsjac.cs$points)))
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[2])]<-"red"
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[3])]<-"green"
points(mdsjac.cs,pch=16,col=colvect)
legend("bottomright",legend = levels(factor(metadata$Factor2)),col=c("blue","red","green"),pch = 16)

plot(mdscyd,type="n",display="sites")
ordipointlabel(mdscyd,display = "sites",add=TRUE)
colvect<-rep("blue",length(rownames(mdscyd$points)))
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[2])]<-"red"
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[3])]<-"green"
points(mdscyd,pch=16,col=colvect)
legend("bottomright",legend = levels(factor(metadata$Factor2)),col=c("blue","red","green"),pch = 16)

plot(mdscyd.cs,type="n",display="sites")
ordipointlabel(mdscyd.cs,display = "sites",add=TRUE)
colvect<-rep("blue",length(rownames(mdscyd.cs$points)))
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[2])]<-"red"
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[3])]<-"green"
points(mdscyd.cs,pch=16,col=colvect)
legend("bottomright",legend = levels(factor(metadata$Factor2)),col=c("blue","red","green"),pch = 16)


stressplot(mdsjac,disjac)
stressplot(mdsjac.cs,disjac.cs)
stressplot(mdscyd,discyd)
stressplot(mdscyd.cs,discyd.cs)


pcoajac <- pcoa(disjac,correction="cailliez")
biplot(pcoajac,t(shared.relab))

pcoajac <- cmdscale(disjac)
plot(pcoajac,type="n",display="sites")
ordipointlabel(pcoajac,display = "sites",add=TRUE)
colvect<-rep("blue",length(rownames(pcoajac$points)))
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[2])]<-"red"
colvect[which(metadata$Factor2==levels(factor(metadata$Factor2))[3])]<-"green"
points(pcoajac,pch=16,col=colvect)
legend("bottomright",legend = levels(factor(metadata$Factor2)),col=c("blue","red","green"),pch = 16)

biplot(pcoajac)

#### endpoint absolute count ####

disjac.ep.abs <- vegdist(t(shared.endpoint.abs),method="jaccard")
discyd.ep.abs <- vegdist(t(shared.endpoint.abs),method="cao") #can handle different sample sizes
mdsjac.ep.abs <- metaMDS(disjac.ep.abs,trace=FALSE,try=50,trymax=2000)
mdscyd.ep.abs <- metaMDS(discyd.ep.abs,trace=FALSE,try=50,trymax=5000)

plot(mdsjac.ep.abs,type="n",display="sites")
ordipointlabel(mdsjac.ep.abs,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdsjac.ep.abs$points)))
colvect[which(metadata.no.inoculum$Factor2==levels(factor(metadata.no.inoculum$Factor2))[2])]<-"red"
pchvect <- rep(15,length(rownames(mdsjac.ep.abs$points)))
pchvect[which(metadata.no.inoculum$Factor1==levels(factor(metadata.no.inoculum$Factor1))[2])]<-16
pchvect[which(metadata.no.inoculum$Factor1==levels(factor(metadata.no.inoculum$Factor1))[3])]<-17
points(mdsjac.ep.abs,pch=pchvect,col=colvect,cex=1.5)
legend("bottomleft",
       legend = paste(rep(levels(factor(metadata.no.inoculum$Factor1)),2),
                      levels(factor(metadata.no.inoculum$Factor2))),
       col=rep(c("blue","red"),3),pch = rep(15:17,2))


plot(mdscyd.ep.abs,type="n",display="sites")
ordipointlabel(mdscyd.ep.abs,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdsjac.ep.abs$points)))
colvect[which(metadata.no.inoculum$Factor2==levels(factor(metadata.no.inoculum$Factor2))[2])]<-"red"
pchvect <- rep(15,length(rownames(mdsjac.ep.abs$points)))
pchvect[which(metadata.no.inoculum$Factor1==levels(factor(metadata.no.inoculum$Factor1))[2])]<-16
pchvect[which(metadata.no.inoculum$Factor1==levels(factor(metadata.no.inoculum$Factor1))[3])]<-17
points(mdscyd.ep.abs,pch=pchvect,col=colvect,cex=1.5)
legend("bottomleft",
       legend = paste(rep(levels(factor(metadata.no.inoculum$Factor1)),2),
                      levels(factor(metadata.no.inoculum$Factor2))),
       col=rep(c("blue","red"),3),pch = rep(15:17,2))



pcoajac <- pcoa(disjac.ep.abs,correction="cailliez")
biplot(pcoajac,t(shared.endpoint.abs))

pcoajac.ep.abs <- cmdscale(disjac.ep.abs)
plot(pcoajac.ep.abs,type="n",display="sites")
ordipointlabel(pcoajac.ep.abs,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdsjac.ep.abs$points)))
colvect[which(metadata.no.inoculum$Factor2==levels(factor(metadata.no.inoculum$Factor2))[2])]<-"red"
pchvect <- rep(15,length(rownames(mdsjac.ep.abs$points)))
pchvect[which(metadata.no.inoculum$Factor1==levels(factor(metadata.no.inoculum$Factor1))[2])]<-16
pchvect[which(metadata.no.inoculum$Factor1==levels(factor(metadata.no.inoculum$Factor1))[3])]<-17
points(pcoajac.ep.abs,pch=pchvect,col=colvect,cex=1.5)
legend("bottom",
       legend = paste(rep(levels(factor(metadata.no.inoculum$Factor1)),2),
                      levels(factor(metadata.no.inoculum$Factor2))),
       col=rep(c("blue","red"),3),pch = rep(15:17,2))

#### BS timeseries ordination ####
shared.bs.abs <- as.data.frame(shared.absolute) %>% 
  dplyr::select(which(grepl(pattern = "BS",colnames(shared.t.ns)))) 
metadata.bs <- metadata[metadata$SampleName %in% colnames(shared.bs.abs),]

disjac.bs.abs <- vegdist(t(shared.bs.abs),method="jaccard")
discyd.bs.abs <- vegdist(t(shared.bs.abs),method="cao") #can handle different sample sizes
mdsjac.bs.abs <- metaMDS(disjac.bs.abs,trace=FALSE,try=50,trymax=2000)
mdscyd.bs.abs <- metaMDS(discyd.bs.abs,trace=FALSE,try=50,trymax=5000)


plot(mdsjac.bs.abs,type="n",display="sites")
ordipointlabel(mdsjac.bs.abs,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdsjac$points)))
colvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[2])]<-"red"
colvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[3])]<-"green"
points(mdsjac.bs.abs,pch=16,col=colvect,cex=1.5)
legend("topright",
       legend = levels(factor(metadata.bs$Factor2)),
       col=c("blue","red","green"),pch = 16)


plot(mdscyd.bs.abs,type="n",display="sites")
ordipointlabel(mdscyd.bs.abs,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdsjac$points)))
colvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[2])]<-"red"
colvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[3])]<-"green"
points(mdscyd.bs.abs,pch=16,col=colvect,cex=1.5)
legend("topright",
       legend = levels(factor(metadata.bs$Factor2)),
       col=c("blue","red","green"),pch = 16)

stressplot(mdsjac.bs.abs,disjac.bs.abs)
stressplot(mdscyd.bs.abs,discyd.bs.abs)


pcoajac.bs.abs <- pcoa(disjac.bs.abs,correction="cailliez")
pcoacyd.bs.abs <- pcoa(discyd.bs.abs,correction="cailliez")
biplot(pcoajac.bs.abs,t(shared.bs.abs))
biplot(pcoacyd.bs.abs,t(shared.bs.abs))

pcoajac.bs.abs.cmd <- cmdscale(disjac.bs.abs,eig=TRUE)
eigvalpcoajac.bs.abs <- scale(abs(pcoajac.bs.abs.cmd$eig),
                              center=FALSE,
                              scale=sum(abs(pcoajac.bs.abs.cmd$eig)))


plot(pcoajac.bs.abs.cmd$points,type="n",
     xlab=paste("PCoA axis 1 (",signif(eigvalpcoajac.bs.abs[1]*100,4),"%)"),
     ylab=paste("PCoA axis 2 (",signif(eigvalpcoajac.bs.abs[2]*100,4),"%)"))
#ordipointlabel(pcoajac.bs.abs.cmd,display="sites",add=TRUE)
colvect <- rep("green",length(rownames(mdsjac$points)))
colvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[2])]<-"blue"
colvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[3])]<-"red"

pchvect <- rep(15,length(rownames(mdscyd.abs.ni$points)))
pchvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[1]&grepl("1",metadata.bs$SampleName))]<-0
pchvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[2]&grepl("2",metadata.bs$SampleName))]<-15
pchvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[2]&grepl("1",metadata.bs$SampleName))]<-0
pchvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[3]&grepl("2",metadata.bs$SampleName))]<-15
pchvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[3]&grepl("1",metadata.bs$SampleName))]<-0

points(pcoajac.bs.abs.cmd$points,
       col=alpha(colvect,0.6),pch=pchvect,cex=1.5)

legend("top",
       legend = c(expression("Inoculum"),
                  expression(paste(BS-NH[4]^"+",", Day44",sep="")),
                  expression(paste("BS-NH"[4]^"+",", Day80")),
                  expression(paste("BS-N"[2],", Day44")),
                  expression(paste("BS-N"[2],", Day80"))),
       col=alpha(c("green","red","red","blue","blue")),pch = c(15,15,0,15,0),
       bty="n")

tiff(filename="BSabscmdjacPcoa_redone.tiff",width =5,height=5,units="in",pointsize = 10,family="sans",res = 300)
plot(pcoajac.bs.abs.cmd$points,type="n",
     xlab=paste("PCoA axis 1 (",signif(eigvalpcoajac.bs.abs[1]*100,4),"%)"),
     ylab=paste("PCoA axis 2 (",signif(eigvalpcoajac.bs.abs[2]*100,4),"%)"))
#ordipointlabel(pcoajac.bs.abs.cmd,display="sites",add=TRUE)
colvect <- rep("green",length(rownames(mdsjac$points)))
colvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[2])]<-"blue"
colvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[3])]<-"red"

pchvect <- rep(15,length(rownames(mdscyd.abs.ni$points)))
pchvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[1]&grepl("1",metadata.bs$SampleName))]<-0
pchvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[2]&grepl("2",metadata.bs$SampleName))]<-15
pchvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[2]&grepl("1",metadata.bs$SampleName))]<-0
pchvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[3]&grepl("2",metadata.bs$SampleName))]<-15
pchvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[3]&grepl("1",metadata.bs$SampleName))]<-0

points(pcoajac.bs.abs.cmd$points,
       col=alpha(colvect,0.6),pch=pchvect,cex=1.5)

legend("top",
       legend = c(expression("Inoculum"),
                  expression(paste(BS-NH[4]^"+",", Day44",sep="")),
                  expression(paste("BS-NH"[4]^"+",", Day80")),
                  expression(paste("BS-N"[2],", Day44")),
                  expression(paste("BS-N"[2],", Day80"))),
       col=alpha(c("green","red","red","blue","blue")),pch = c(15,15,0,15,0),
       bty="n")

dev.off()

biplot(pcoajac)

pcoacyd.bs.abs.cmd <- cmdscale(discyd.bs.abs)
plot(pcoacyd.bs.abs.cmd,type="n",display="sites")
ordipointlabel(pcoacyd.bs.abs.cmd,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdsjac$points)))
colvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[2])]<-"red"
colvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[3])]<-"green"
points(pcoacyd.bs.abs.cmd,pch=16,col=colvect,cex=1.5)
legend("bottom",
       legend = levels(factor(metadata.bs$Factor2)),
       col=c("blue","red","green"),pch = 16)

#### HE timeseries ordination ####
shared.he.abs <- as.data.frame(shared.absolute) %>% 
  dplyr::select(which(grepl(pattern = "HE",colnames(shared.t.ns)))) 
metadata.he <- metadata[metadata$SampleName %in% colnames(shared.he.abs),]

disjac.he.abs <- vegdist(t(shared.he.abs),method="jaccard")
discyd.he.abs <- vegdist(t(shared.he.abs),method="cao") #can handle different sample sizes
mdsjac.he.abs <- metaMDS(disjac.he.abs,trace=FALSE,try=50,trymax=2000)
mdscyd.he.abs <- metaMDS(discyd.he.abs,trace=FALSE,try=50,trymax=5000)


plot(mdsjac.he.abs,type="n",display="sites")
ordipointlabel(mdsjac.he.abs,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdsjac$points)))
colvect[which(metadata.he$Factor2==levels(factor(metadata.he$Factor2))[2])]<-"red"
colvect[which(metadata.he$Factor2==levels(factor(metadata.he$Factor2))[3])]<-"green"
points(mdsjac.he.abs,pch=16,col=colvect,cex=1.5)
legend("topright",
       legend = levels(factor(metadata.he$Factor2)),
       col=c("blue","red","green"),pch = 16)


plot(mdscyd.he.abs,type="n",display="sites")
ordipointlabel(mdscyd.he.abs,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdsjac$points)))
colvect[which(metadata.he$Factor2==levels(factor(metadata.he$Factor2))[2])]<-"red"
colvect[which(metadata.he$Factor2==levels(factor(metadata.he$Factor2))[3])]<-"green"
points(mdscyd.he.abs,pch=16,col=colvect,cex=1.5)
legend("topright",
       legend = levels(factor(metadata.bs$Factor2)),
       col=c("blue","red","green"),pch = 16)

stressplot(mdsjac.he.abs,disjac.he.abs)
stressplot(mdscyd.he.abs,discyd.he.abs)


pcoajac.he.abs <- pcoa(disjac.he.abs,correction="cailliez")
pcoacyd.he.abs <- pcoa(discyd.he.abs,correction="cailliez")
biplot(pcoajac.he.abs,t(shared.he.abs))
biplot(pcoacyd.he.abs,t(shared.he.abs))

pcoajac.he.abs.cmd <- cmdscale(disjac.he.abs)
plot(pcoajac.he.abs.cmd,type="n",display="sites")
ordipointlabel(pcoajac.he.abs.cmd,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdsjac$points)))
colvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[2])]<-"red"
colvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[3])]<-"green"
points(pcoajac.he.abs.cmd,pch=16,col=colvect,cex=1.5)
legend("top",
       legend = levels(factor(metadata.bs$Factor2)),
       col=c("blue","red","green"),pch = 16)

biplot(pcoajac)

pcoacyd.he.abs.cmd <- cmdscale(discyd.he.abs)
plot(pcoacyd.he.abs.cmd,type="n",display="sites")
ordipointlabel(pcoacyd.he.abs.cmd,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdsjac$points)))
colvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[2])]<-"red"
colvect[which(metadata.bs$Factor2==levels(factor(metadata.bs$Factor2))[3])]<-"green"
points(pcoacyd.he.abs.cmd,pch=16,col=colvect,cex=1.5)
legend("bottom",
       legend = levels(factor(metadata.bs$Factor2)),
       col=c("blue","red","green"),pch = 16)

pcoajac.he.abs.cmd <- cmdscale(disjac.he.abs,eig=TRUE)
eigvalpcoajac.he.abs <- scale(abs(pcoajac.he.abs.cmd$eig),
                              center=FALSE,
                              scale=sum(abs(pcoajac.he.abs.cmd$eig)))


plot(pcoajac.he.abs.cmd$points,type="n",
     xlab=paste("PCoA axis 1 (",signif(eigvalpcoajac.he.abs[1]*100,4),"%)"),
     ylab=paste("PCoA axis 2 (",signif(eigvalpcoajac.he.abs[2]*100,4),"%)"))
#ordipointlabel(pcoajac.he.abs.cmd,display="sites",add=TRUE)
colvect <- rep("green",length(rownames(mdsjac$points)))
colvect[which(metadata.he$Factor2==levels(factor(metadata.he$Factor2))[2])]<-"blue"
colvect[which(metadata.he$Factor2==levels(factor(metadata.he$Factor2))[3])]<-"red"

pchvect <- rep(16,length(rownames(mdscyd.abs.ni$points)))
pchvect[which(metadata.he$Factor2==levels(factor(metadata.he$Factor2))[1]&grepl("1",metadata.he$SampleName))]<-1
pchvect[which(metadata.he$Factor2==levels(factor(metadata.he$Factor2))[2]&grepl("2",metadata.he$SampleName))]<-16
pchvect[which(metadata.he$Factor2==levels(factor(metadata.he$Factor2))[2]&grepl("1",metadata.he$SampleName))]<-1
pchvect[which(metadata.he$Factor2==levels(factor(metadata.he$Factor2))[3]&grepl("2",metadata.he$SampleName))]<-16
pchvect[which(metadata.he$Factor2==levels(factor(metadata.he$Factor2))[3]&grepl("1",metadata.he$SampleName))]<-1

points(pcoajac.he.abs.cmd$points,
       col=alpha(colvect,0.6),pch=pchvect,cex=1.5)

legend("top",
       legend = c(expression("Inoculum"),
                  expression(paste(HE-NH[4]^"+",", Day44",sep="")),
                  expression(paste("HE-NH"[4]^"+",", Day80")),
                  expression(paste("HE-N"[2],", Day44")),
                  expression(paste("HE-N"[2],", Day80"))),
       col=alpha(c("green","red","red","blue","blue")),pch = c(16,16,1,16,1),
       bty="n")

tiff(filename="HEabscmdjacPcoa_redone.tiff",width =5,height=5,units="in",pointsize = 10,family="sans",res = 300)
plot(pcoajac.he.abs.cmd$points,type="n",
     xlab=paste("PCoA axis 1 (",signif(eigvalpcoajac.he.abs[1]*100,4),"%)"),
     ylab=paste("PCoA axis 2 (",signif(eigvalpcoajac.he.abs[2]*100,4),"%)"))
#ordipointlabel(pcoajac.he.abs.cmd,display="sites",add=TRUE)
colvect <- rep("green",length(rownames(mdsjac$points)))
colvect[which(metadata.he$Factor2==levels(factor(metadata.he$Factor2))[2])]<-"blue"
colvect[which(metadata.he$Factor2==levels(factor(metadata.he$Factor2))[3])]<-"red"

pchvect <- rep(16,length(rownames(mdscyd.abs.ni$points)))
pchvect[which(metadata.he$Factor2==levels(factor(metadata.he$Factor2))[1]&grepl("1",metadata.he$SampleName))]<-1
pchvect[which(metadata.he$Factor2==levels(factor(metadata.he$Factor2))[2]&grepl("2",metadata.he$SampleName))]<-16
pchvect[which(metadata.he$Factor2==levels(factor(metadata.he$Factor2))[2]&grepl("1",metadata.he$SampleName))]<-1
pchvect[which(metadata.he$Factor2==levels(factor(metadata.he$Factor2))[3]&grepl("2",metadata.he$SampleName))]<-16
pchvect[which(metadata.he$Factor2==levels(factor(metadata.he$Factor2))[3]&grepl("1",metadata.he$SampleName))]<-1

points(pcoajac.he.abs.cmd$points,
       col=alpha(colvect,0.6),pch=pchvect,cex=1.5)

legend("top",
       legend = c(expression("Inoculum"),
                  expression(paste(HE-NH[4]^"+",", Day44",sep="")),
                  expression(paste("HE-NH"[4]^"+",", Day80")),
                  expression(paste("HE-N"[2],", Day44")),
                  expression(paste("HE-N"[2],", Day80"))),
       col=alpha(c("green","red","red","blue","blue")),pch = c(16,16,1,16,1),
       bty="n")
dev.off()



#### SS timeseries ordination ####
shared.ss.abs <- as.data.frame(shared.absolute) %>% 
  dplyr::select(which(grepl(pattern = "SS",colnames(shared.t.ns)))) 
metadata.ss <- metadata[metadata$SampleName %in% colnames(shared.ss.abs),]

disjac.ss.abs <- vegdist(t(shared.ss.abs),method="jaccard")
discyd.ss.abs <- vegdist(t(shared.ss.abs),method="cao") #can handle different sample sizes
mdsjac.ss.abs <- metaMDS(disjac.ss.abs,trace=FALSE,try=50,trymax=2000)
mdscyd.ss.abs <- metaMDS(discyd.ss.abs,trace=FALSE,try=50,trymax=5000)


plot(mdsjac.ss.abs,type="n",display="sites")
ordipointlabel(mdsjac.ss.abs,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdsjac$points)))
colvect[which(metadata.ss$Factor2==levels(factor(metadata.ss$Factor2))[2])]<-"red"
points(mdsjac.ss.abs,pch=16,col=colvect,cex=1.5)
legend("topright",
       legend = levels(factor(metadata.ss$Factor2)),
       col=c("blue","red"),pch = 16)


plot(mdscyd.ss.abs,type="n",display="sites")
ordipointlabel(mdscyd.ss.abs,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdsjac$points)))
colvect[which(metadata.ss$Factor2==levels(factor(metadata.ss$Factor2))[2])]<-"red"
points(mdscyd.ss.abs,pch=16,col=colvect,cex=1.5)
legend("topright",
       legend = levels(factor(metadata.ss$Factor2)),
       col=c("blue","red"),pch = 16)

stressplot(mdsjac.ss.abs,disjac.ss.abs)
stressplot(mdscyd.ss.abs,discyd.ss.abs)


pcoajac.ss.abs <- pcoa(disjac.ss.abs,correction="cailliez")
pcoacyd.ss.abs <- pcoa(discyd.ss.abs,correction="cailliez")
biplot(pcoajac.ss.abs,t(shared.ss.abs))
biplot(pcoacyd.ss.abs,t(shared.ss.abs))

pcoajac.ss.abs.cmd <- cmdscale(disjac.ss.abs)
plot(pcoajac.ss.abs.cmd,type="n",display="sites")
ordipointlabel(pcoajac.ss.abs.cmd,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdsjac$points)))
colvect[which(metadata.ss$Factor2==levels(factor(metadata.ss$Factor2))[2])]<-"red"
points(pcoajac.ss.abs.cmd,pch=16,col=colvect,cex=1.5)
legend("top",
       legend = levels(factor(metadata.ss$Factor2)),
       col=c("blue","red"),pch = 16)

biplot(pcoajac.ss.abs.cmd)

pcoacyd.ss.abs.cmd <- cmdscale(discyd.ss.abs)
plot(pcoacyd.ss.abs.cmd,type="n",display="sites")
ordipointlabel(pcoacyd.ss.abs.cmd,display = "sites",add=TRUE)
colvect <- rep("blue",length(rownames(mdsjac$points)))
colvect[which(metadata.ss$Factor2==levels(factor(metadata.ss$Factor2))[2])]<-"red"
points(pcoacyd.ss.abs.cmd,pch=16,col=colvect,cex=1.5)
legend("bottom",
       legend = levels(factor(metadata.ss$Factor2)),
       col=c("blue","red"),pch = 16)


pcoajac.ss.abs.cmd <- cmdscale(disjac.ss.abs,eig=TRUE)
eigvalpcoajac.ss.abs <- scale(abs(pcoajac.ss.abs.cmd$eig),
                              center=FALSE,
                              scale=sum(abs(pcoajac.ss.abs.cmd$eig)))


plot(pcoajac.ss.abs.cmd$points,type="n",
     xlab=paste("PCoA axis 1 (",signif(eigvalpcoajac.ss.abs[1]*100,4),"%)"),
     ylab=paste("PCoA axis 2 (",signif(eigvalpcoajac.ss.abs[2]*100,4),"%)"))
#ordipointlabel(pcoajac.ss.abs.cmd,display="sites",add=TRUE)
colvect <- rep("blue",length(rownames(pcoajac.ss.abs.cmd$points)))
colvect[which(metadata.ss$Factor2==levels(factor(metadata.ss$Factor2))[2])]<-"red"


pchvect <- rep(17,length(rownames(pcoajac.ss.abs.cmd$points)))
pchvect[which(metadata.ss$Factor2==levels(factor(metadata.ss$Factor2))[1]&grepl("1",metadata.ss$SampleName))]<-2
pchvect[which(metadata.ss$Factor2==levels(factor(metadata.ss$Factor2))[2]&grepl("2",metadata.ss$SampleName))]<-17
pchvect[which(metadata.ss$Factor2==levels(factor(metadata.ss$Factor2))[2]&grepl("1",metadata.ss$SampleName))]<-2


points(pcoajac.ss.abs.cmd$points,
       col=alpha(colvect,0.6),pch=pchvect,cex=1.5)

legend("top",
       legend = c(expression(paste(SS-NH[4]^"+",", Day53",sep="")),
                  expression(paste("SS-NH"[4]^"+",", Day89")),
                  expression(paste("SS-N"[2],", Day53")),
                  expression(paste("SS-N"[2],", Day89"))),
       col=alpha(c("red","red","blue","blue"),0.6),pch = c(17,2,17,2),
       bty="n")

tiff(filename="SSabscmdjacPcoa_redone.tiff",width =5,height=5,units="in",pointsize = 10,family="sans",res = 300)
plot(pcoajac.ss.abs.cmd$points,type="n",
     xlab=paste("PCoA axis 1 (",signif(eigvalpcoajac.ss.abs[1]*100,4),"%)"),
     ylab=paste("PCoA axis 2 (",signif(eigvalpcoajac.ss.abs[2]*100,4),"%)"))
#ordipointlabel(pcoajac.ss.abs.cmd,display="sites",add=TRUE)
colvect <- rep("blue",length(rownames(pcoajac.ss.abs.cmd$points)))
colvect[which(metadata.ss$Factor2==levels(factor(metadata.ss$Factor2))[2])]<-"red"


pchvect <- rep(17,length(rownames(pcoajac.ss.abs.cmd$points)))
pchvect[which(metadata.ss$Factor2==levels(factor(metadata.ss$Factor2))[1]&grepl("1",metadata.ss$SampleName))]<-2
pchvect[which(metadata.ss$Factor2==levels(factor(metadata.ss$Factor2))[2]&grepl("2",metadata.ss$SampleName))]<-17
pchvect[which(metadata.ss$Factor2==levels(factor(metadata.ss$Factor2))[2]&grepl("1",metadata.ss$SampleName))]<-2


points(pcoajac.ss.abs.cmd$points,
       col=alpha(colvect,0.6),pch=pchvect,cex=1.5)

legend("top",
       legend = c(expression(paste(SS-NH[4]^"+",", Day53",sep="")),
                  expression(paste("SS-NH"[4]^"+",", Day89")),
                  expression(paste("SS-N"[2],", Day53")),
                  expression(paste("SS-N"[2],", Day89"))),
       col=alpha(c("red","red","blue","blue"),0.6),pch = c(17,2,17,2),
       bty="n")

dev.off()