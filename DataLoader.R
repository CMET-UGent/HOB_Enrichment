#### load required packages #### 
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
#### custom functions ####
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
    tax.final$Species <- paste(tax.final$Genus,tax.final$Species) #just the species is not very sensible
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



#extrafont::font_import()

#### load the data ####
#read shared file (OTU table)
fldr <- "NGS881_Xiaona"
fn <- "stability"
shared <- fread(input = "NGS881_Xiaona/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared",
                header=TRUE) #replace here with the filename and location of your mothur shared file
shared <- as.data.frame(shared) # fread reads in as tibble, convert to data.frame for downstream processing
desgroups <- shared$Group # select sample names
shared.x <- shared[,4:ncol(shared)] # remove all non-count data from the OTU table
rownames(shared.x) <- desgroups
shared.t <- as.data.frame(t(shared.x))
shared.t.ns <- shared.t[which(rowSums(shared.t)!=1),]
sharedns <- data.frame(label=rep(0.03,ncol(shared.t.ns)),
                       Group=colnames(shared.t.ns),
                       numOtus=nrow(shared.t.ns),t(shared.t.ns))
# read taxonomy file
otutaxonomy <- fread(input = "NGS881_Xiaona/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy",header=TRUE) # replace here with the filename and location of your mothur taxonomy file

taxonomy.spl <- preformattax(otutaxonomy)
taxonomy.np <- taxonomy.spl %>% dplyr::select(-dplyr::contains("Prob"))

#filter taxonomy
taxonomy.np.ns <- taxonomy.np[which(rownames(taxonomy.np) 
                                    %in% rownames(shared.t.ns)),]
metadata.tibble <- readxl::read_excel("NGS881_Xiaona/MetaData_5.xlsx",sheet="ForR") #set location of your own metadata.xlsx
factdescs <- readxl::read_excel("NGS881_Xiaona/MetaData_5.xlsx",sheet="FactDesc")
metadata <- as.data.frame(metadata.tibble) #to avoid warnings/errors with rownames
if(is.numeric(metadata$SampleName)) #check for fully numeric sample names
{
  metadata$SampleName <- paste(dataname,metadata$SampleName,sep="")
}
rownames(metadata) <- metadata$SampleName

colnames(shared.t.ns) <- plyr::mapvalues(colnames(shared.t.ns),
                                         from=as.character(metadata$Code),
                                         to=as.character(metadata$SampleName)) #rename

shared.t.ns <- shared.t.ns %>% select(-c(`SS0-a`,`SS0-b`,`SS0-c`))
metadata <- metadata %>%  filter(!SampleName%in%c("SS0-a","SS0-b","SS0-c"))
rownames(metadata) <- metadata$SampleName
metadata.smpdat <- sample_data(metadata)

# read celcount info
cellconcentations <- read.csv("Cellcounts.csv",header = FALSE, col.names = c("Sample","ConcentrationpermL"))

#### phyloseq constructors ####
otumat.ns <- as.matrix(shared.t.ns)
taxmat.ns <- as.matrix(taxonomy.np.ns)
OTU       <- otu_table(otumat.ns,taxa_are_rows = TRUE)
TAX       <- tax_table(taxmat.ns)
physeqobj <- phyloseq(OTU,TAX)
physeqobj.meta <- merge_phyloseq(physeqobj,metadata.smpdat)

#### filtering according to WNWN ####
## prevalence = (fraction of samples in which an OTU is observed minimum 1 time)
minobs=1
prevalence <- apply(as.matrix(shared.t.ns),1,function(x,minobs){sum(x>=minobs)},minobs)/ncol(shared.t.ns)
prevalencefilter <- prevalence>0.05
sharedminsingletonwnwn <- shared.t.ns[prevalencefilter,]
##Read counts should exceed 0.5 times the number of samples
sharedfilteredwnwn <- sharedminsingletonwnwn[rowSums(sharedminsingletonwnwn)>0.5*ncol(sharedminsingletonwnwn),]
# deseq normalise ==> very dependent upon design
metadata$factor1 <- factor(metadata$Factor1)
metadata$factor2 <- factor(metadata$Factor2)
deseqdata <- as.matrix(sharedfilteredwnwn +1)

#### heatmaps without standardization ####
relabsphyseq <- transform_sample_counts(physeqobj, function(x) x/sum(x)) # create relative abudances
Top24OTUs.relab <- names(sort(taxa_sums(relabsphyseq), TRUE)[1:24])      # select top 24 OTUs in relative abundance
physeqobj24.relab  <- prune_taxa(Top24OTUs.relab, relabsphyseq)          # extract the top 24 OTUs
relab24mat <- as.matrix(otu_table(physeqobj24.relab))                    # Create a matrix with the selected relative abundances
rownames(relab24mat) <- as.character(as.data.frame(tax_table(physeqobj24.relab)@.Data)$Genus)
my_palette <- colorRampPalette(c("black","green" ,"red"))(n = 409)
# this colors can be selected from a hexadecimal color picker
# e.g. https://www.google.be/search?hl=en&q=hex+color+picker&gws_rd=ssl

#heatmap.2(relab24mat,col=brewer.pal(9,name="Reds"))



topn <- 20
wnwnsh <- sharedfilteredwnwn
wnwnsh.prop <- decostand(wnwnsh,method = "total",MARGIN=2)
shns <- shared.t.ns
shns.prop <- decostand(shns,method = "total",MARGIN=2)
shns.log <- decostand(shns,method = "log",logbase=10)
shns.cs <- shns.prop*as.numeric(quantile(colSums(shared.t.ns),.75)) #other quantile then max see: #boxplot(colSums(shared.t.ns))
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
wnwnsh.top.clust <- hclust(vegdist(t(wnwnsh.top),"jaccard",binary = FALSE),method="ward.D2")

shared.relab <- decostand(shared.t.ns,method="total",MARGIN=2)
shared.cs <- shared.relab*min(colSums(shared.t.ns)) #common scale strategy

shared.absolute <- round(t(t(shared.relab)*cellconcentations$ConcentrationpermL))

shared.no.inoculum <- shared.t.ns %>% dplyr::select(-contains("0"))
shared.no.inoculum.abs <- as.data.frame(shared.absolute) %>% dplyr::select(-contains("0"))
metadata.no.inoculum <- metadata %>% dplyr::filter(!grepl("0",rownames(metadata)))
rownames(metadata.no.inoculum) <- rownames(metadata)[!grepl("0",rownames(metadata))]