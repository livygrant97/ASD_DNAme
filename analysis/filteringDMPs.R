################################################################
# Characterising Autosomal sex differences between males and females
# Olivia Grant
###############################################################

#im loading a lot of packages here so just using l apply
packages <- c("missMethyl", "minfi", "cluster","IlluminaHumanMethylationEPICanno.ilm10b2.hg19","dplyr",
"bigmelon","clusterProfiler","DMRcate","ggplot2","factoextra","tidyverse","IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
lapply(packages, library, character.only = TRUE)


###############################################################
#### loading and cleaning the discovery data ####
###############################################################
dim(dmp)

data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
annotationTable <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
annotationTable<-as.data.frame(annotationTable)
#get probes on x and y chromosome
sexProbes<-as.character(annotationTable$Name[annotationTable$chr %in% c("chrX","chrY")])
SNPAnnotations<-read.delim("/home/og16379/diff_cpg_fm/data/EPIC.hg19.manifest.tsv")
SNPAssociatedCpG<- SNPAnnotations %>% dplyr::filter(MASK_general==TRUE)
chProbes <- read.table("/home/og16379/diff_cpg_fm/data/chProbes.txt")

#performing FDR test
dmp$FDR<-p.adjust(dmp$P.Value,method="fdr")
#subset those probes which FDR <0.05
Passed_FDR<- subset(dmp,FDR < 0.05)
#load annotations
#remove probes not on autosomes
sigCpG_Autosomes<- Passed_FDR[!Passed_FDR$Row.names %in% sexProbes,]
#read in annotations for SNP related probes remove them from analysis
dasenAutosomeCleaned<- sigCpG_Autosomes[!sigCpG_Autosomes$Row.names %in%  SNPAssociatedCpG$probeID,]
# downloaded the list of cross hybridizing probes from here https://pubmed.ncbi.nlm.nih.gov/23314698/
dasenAutosomeCleaned<-dasenAutosomeCleaned[!(dasenAutosomeCleaned$Row.names
  %in% chProbes$V1),]
paste(formatC(nrow(SNPAssociatedCpG),big.mark = ','),'Number of SNP associated probes removed')
paste(formatC(nrow(chProbes),big.mark = ','),'Number of cross hybridizing probes removed')
save(dasenAutosomeCleaned,file="/home/og16379/diff_cpg_fm/data/dasenAutosomeCleanedFirstRound.Rata")
dim(dasenAutosomeCleaned)


###############################################################
####  Filtering list to get strongest differences with discovery ####
###############################################################

# then we wanted to filter this list down to select the CpGs which showed the strongest differences
# between males and females
# first, detirmine what sex methylation is higher in for each probe using the average beta values
# add a column called $Higher_Met_In which indicates which sex methylation is higher in
dasenAutosomeCleaned$Higher_Met_In<-"Male"
dasenAutosomeCleaned$Higher_Met_In [dasenAutosomeCleaned$F_AVG >dasenAutosomeCleaned$M_AVG] <- "Female"
#split for sex so you have all probes regardless of deltaBeta. incase want to change threshold
splitCpG_sig<-split(dasenAutosomeCleaned,dasenAutosomeCleaned$Higher_Met_In)
femalesSigCpG<-splitCpG_sig$Female
malesSigCpG<-splitCpG_sig$Male
# we used threshold of delta beta diff greater than 0.05
# make object for each for further analysis
sigProbesFemales<- subset(femalesSigCpG, deltaBeta > 0.05  | deltaBeta <  -0.05)
sigProbesMales<- subset(malesSigCpG, deltaBeta > 0.05 | deltaBeta < -0.05)
allSigProbes<-rbind(sigProbesMales,sigProbesFemales)
save(sigProbesFemales,file="/home/og16379/diff_cpg_fm/data/sigProbesFemales.RData")
save(sigProbesMales,file="/home/og16379/diff_cpg_fm/data/sigProbesMales.RData")
save(allSigProbes,file="/home/og16379/diff_cpg_fm/data/allSigProbes.RData")


###############################################################
#### loading and cleaning the validation data ####
###############################################################

dim(dmp2)

#performing FDR test
dmp2$FDR<-p.adjust(dmp2$P.Value,method="fdr")
#subset those probes which FDR <0.05
Passed_FDR2<- subset(dmp2,FDR < 0.05)
#load annotations
#remove probes not on autosomes
sigCpG_Autosomes2<- Passed_FDR2[!Passed_FDR2$Row.names %in% sexProbes,]
#read in annotations for SNP related probes remove them from analysis
dasenAutosomeCleaned2<- sigCpG_Autosomes2[!sigCpG_Autosomes2$Row.names %in%  SNPAssociatedCpG$probeID,]
# downloaded the list of cross hybridizing probes from here https://pubmed.ncbi.nlm.nih.gov/23314698/
dasenAutosomeCleaned2<-dasenAutosomeCleaned2[!(dasenAutosomeCleaned2$Row.names
  %in% chProbes$V1),]
paste(formatC(nrow(SNPAssociatedCpG),big.mark = ','),'Number of SNP associated probes removed')
paste(formatC(nrow(chProbes),big.mark = ','),'Number of cross hybridizing probes removed')
save(dasenAutosomeCleaned2,file="/home/og16379/diff_cpg_fm/data/dasenAutosomeCleanedSecondRound.Rata")
dim(dasenAutosomeCleaned2)

###############################################################
####  Filtering list to get strongest differences with validation ####
###############################################################

# then we wanted to filter this list down to select the CpGs which showed the strongest differences
# between males and females
# first, detirmine what sex methylation is higher in for each probe using the average beta values
# add a column called $Higher_Met_In which indicates which sex methylation is higher in

dasenAutosomeCleaned2$Higher_Met_In<-"Male"
dasenAutosomeCleaned2$Higher_Met_In [dasenAutosomeCleaned2$F_AVG >dasenAutosomeCleaned2$M_AVG] <- "Female"
#split for sex so you have all probes regardless of deltaBeta. incase want to change threshold
splitCpG_sig2<-split(dasenAutosomeCleaned2,dasenAutosomeCleaned2$Higher_Met_In)
femalesSigCpG2<-splitCpG_sig2$Female
malesSigCpG2<-splitCpG_sig2$Male
# we used threshold of delta beta diff greater than 0.05
# make object for each for further analysis
sigProbesFemales2<- subset(femalesSigCpG2, deltaBeta > 0.05  | deltaBeta < -0.05)
sigProbesMales2<- subset(malesSigCpG2, deltaBeta > 0.05 | deltaBeta < -0.05)
allSigProbes2<-rbind(sigProbesMales2,sigProbesFemales2)
save(sigProbesFemales2,file="/home/og16379/diff_cpg_fm/data/sigProbesFemales2.RData")
save(sigProbesMales2,file="/home/og16379/diff_cpg_fm/data/sigProbesMales2.RData")
save(allSigProbes2,file="/home/og16379/diff_cpg_fm/data/allSigProbes2.RData")


###########################################################################
####  now lets see whats replicated in the larger list and do pathway analysis!!!! ####
###########################################################################
# which results consistent across data sets
replicatedCleaned <- which(dasenAutosomeCleaned$Row.names %in% dasenAutosomeCleaned2$Row.names)
replicatedCleaned <- dasenAutosomeCleaned2[replicatedCleaned,]

# character vector of all CpG sites tested.
all.cpg <- as.character(rownames(annotationTable))
# character vector of significant CpG sites to test for GO term enrichment
sig.cpg<-as.character(replicatedCleaned$Row.names))

# perform gene ontology testing collection GO
sexAssociatedSitesPathwaysExtended <- gometh(sig.cpg,all.cpg,collection="GO",array.type="EPIC")
#get signficiant pathways
table(sexAssociatedSitesPathwaysExtended$FDR<0.05)
#table of top GO results
topGSA(sexAssociatedSitesPathwaysExtended)
# 20 = TRUE
sigPathwaysExtendedGO <- topGSA(sexAssociatedSitesPathwaysExtended,n=20)
write.csv(sigPathwaysExtendedGO,file="/home/og16379/diff_cpg_fm/data/sexAssociatedSitesPathwaysExtendedGO.csv")

#repeat for KEGG collection
sexAssociatedSitesPathwaysExtendedKEGG <- gometh(sig.cpg,all.cpg,collection="KEGG",array.type="EPIC")
#get signficiant pathways
table(sexAssociatedSitesPathwaysExtendedKEGG$FDR<0.05)
# 20 = TRUE
#table of top GO results
sigPathwaysExtendedKEGG <- topGSA(sexAssociatedSitesPathwaysExtendedKEGG,n=20)
write.csv(sigPathwaysExtendedKEGG,file="/home/og16379/diff_cpg_fm/data/sexAssociatedSitesPathwaysExtendedKegg.csv")


###########################################################################
####  now lets see whats replicated in the decreased list  ####
###########################################################################

replicatedResults <- which(allSigProbes2$Row.names %in% allSigProbes$Row.names)
replicatedResults<-allSigProbes[replicatedResults,]
dim(replicatedResults)
#[1] 409  25
write.csv(replicatedResults, file="/home/og16379/diff_cpg_fm/data/sadmpReplicatedAdditionalFile.csv")
save(replicatedResults,file="/home/og16379/diff_cpg_fm/data/replicatedResults.RData")

replicatedDiscovery <- subset(allSigProbes,Row.names %in% replicatedResults$Row.names)
replicatedValidation<- subset(allSigProbes2,Row.names %in% replicatedResults$Row.names)
head(sign(replicatedDiscovery$deltaBeta),n=500)
head(sign(replicatedValidation$deltaBeta),n=500)

 #directions of effect are all the same
replicatedResultsSplit<-split(replicatedResults,replicatedResults$Higher_Met_In)
repFemalesSigCpG<-replicatedResultsSplit$Female
repMalesSigCpG<-replicatedResultsSplit$Male
save(repFemalesSigCpG,file="/home/og16379/diff_cpg_fm/data/repFemalesSigCpG.RData")
save(repMalesSigCpG,file="/home/og16379/diff_cpg_fm/data/repMalesSigCpG.RData")

###########################################################################
####   do pathway analysis ####
###########################################################################

# character vector of all CpG sites tested.
# all.cpg <- as.character(rownames(annotationTable))
# character vector of significant CpG sites to test for GO term enrichment
sig.cpg<-as.character(replicatedResults$Row.names)

# perform gene ontology testing collection GO
sexAssociatedSitesPathwaysDecreased<- gometh(sig.cpg,all.cpg,collection="GO",array.type="EPIC")
#get signficiant pathways
table(sexAssociatedSitesPathwaysDecreased$FDR<0.05)
#table of top GO results
topGSA(sexAssociatedSitesPathwaysDecreased)
# all false

#repeat for KEGG collection
sexAssociatedSitesPathwaysDecreasedKEGG <- gometh(sig.cpg,all.cpg,collection="KEGG",array.type="EPIC")
#get signficiant pathways
table(sexAssociatedSitesPathwaysDecreasedKEGG$FDR<0.05)
# all false


#################################################################################
#### set colours now for plotting the data ####
#################################################################################

#colours we use for plotting throughout the manuscript
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
maleColour<-cbbPalette[7]
femaleColour<-cbbPalette[6]
