################################################################
# Characterising Autosomal sex differences between males and females
# Olivia Grant
###############################################################

#below is the code to perform the tf motif enrichment

################################################################################
# libraries
################################################################################
if(!require("MotifDb", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("MotifDb")
}
library(MotifDb)


if(!require("PWMEnrich", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("PWMEnrich")
}
library(PWMEnrich)



if(!require("PWMEnrich.Hsapiens.background", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("PWMEnrich.Hsapiens.background")
}
library(PWMEnrich.Hsapiens.background)


if(!require("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
}
library(BSgenome.Hsapiens.UCSC.hg19)


#make granges for sex associated probes
gfile <- openfn.gds('/storage/st05d/Exeter/UnderSocMeth/USM1/USM_WF.gds',readonly=T)
##re code the gds2mset function from big melon with right annotation
##re code the gds2mset function from big melon with right annotation
gds2mset <- function(gds, i, j, anno = NULL){
     x <- gds
     if(!is.null(anno)){
         if(!anno %in% c("27k", "450k", "epic")){
         stop("anno needs to be: \'27k\', \'450k\', \'epic\'")
         }
     }
     M <- x[i = i, j = j,   "methylated", name = TRUE, drop = FALSE]
     U <- x[i = i, j = j, "unmethylated", name = TRUE, drop = FALSE]
     pd <- pData(x)[j, , drop = FALSE]
     rownames(pd) <- colnames(x)[j]
     #    pd <- annotatedDataFrameFrom(object = as.matrix(pd), byrow = TRUE)
     if(!is.null(anno)){
         if(anno == "27k"){
             anno <- c("IlluminaHumanMethylation27k", "ilmn12.hg19")
         } else if(anno == "450k"){
             anno <- c("IlluminaHumanMethylation450k", "ilmn12.hg19")
         } else if(anno == "epic"){
             anno <- c("IlluminaHumanMethylationEPIC", "ilm10b4.hg19")
         } else if(anno == "unknown"){
             anno <- c("Unknown", "Unknown")
         }
     }
     # Guess Array Type - will not get correct array if performed on subset.
     if(is.null(anno)){
         nr <- nrow(fData(x))
         # Will guess array type based on number of rows, will fail on subsets!
         if(nr > 50000 & nr < 500000){
             anno <- c("IlluminaHumanMethylation450k", "ilmn12.hg19")
         } else if(nr >= 500000){
             anno <- c("IlluminaHumanMethylationEPIC", "ilm10b4.hg19")
         } else if(nr <=50000){
             anno <- c("IlluminaHumanMethylation27k", "ilmn12.hg19")
         }
     }
     names(anno) <- c("array", "annotation")
     out <- MethylSet(Meth = M, Unmeth = U, colData = pd, annotation = anno)
     out@preprocessMethod <- c(
         rg.norm = "Converted from gdsfmt to MethylSet (bigmelon)",
         minfi = as.character(packageVersion("minfi")),
         manifest = NA #packageVersion(getManifest(anno))
         )
     out
 }
#convert to mset
gfile<-gds2mset(gfile)
#map it to the genome
msetMapped<-mapToGenome(gfile)
#convert that to a g range object
msetMapped<-granges(msetMapped)

#subset that by saDMPs
sigProbeMapped<-subset(msetMapped,msetMapped@ranges@NAMES %in% replicatedResults$Row.names)
rtracklayer::export.bed(sigProbeMapped,"~/res.bed")

femaleProbeMapped<-subset(msetMapped,msetMapped@ranges@NAMES %in% repFemalesSigCpG$Row.names)
maleProbeMapped<-subset(msetMapped,msetMapped@ranges@NAMES %in% repMalesSigCpG$Row.names)
rtracklayer::export.bed(femaleProbeMapped,"~/resFemaleRep.bed")
rtracklayer::export.bed(maleProbeMapped,"~/resMaleRep.bed")



################################################################################
# get data
################################################################################
# load in data
data(PWMEnrich.Hsapiens.background)
data(PWMLogn.hg19.MotifDb.Hsap)
# register cores
registerCoresPWMEnrich(30)

# extend the GRanges by 100bp
femaleProbeMapped<-resize(femaleProbeMapped,fix="start",width=width(femaleProbeMapped)+100)
maleProbeMapped<-resize(maleProbeMapped,fix="start",width=width(maleProbeMapped)+100)

# make dna string sets
femaleDMPseq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, femaleGRanges)
maleDMPseq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, maleGRanges)

# get PWMS
hg19_PWMs <- PWMLogn.hg19.MotifDb.Hsap$pwms


################################################################################
# build backgrounds
################################################################################
if(file.exists("/home/og16379/diff_cpg_fm/data/male_background.RData")){
  load("/home/og16379/diff_cpg_fm/data/male_background.RData")
} else{
  male_seq_bg <- makeBackground(motifs = PWMLogn.hg19.MotifDb.Hsap$pwms, bg.seq=maleDMPseq,
                                              type="logn", algorithm="human",  verbose=TRUE)

  save(male_seq_bg, file="/home/og16379/diff_cpg_fm/data/male_background.RData")
}


if(file.exists("/home/og16379/diff_cpg_fm/data/female_background.RData")){
  load("/home/og16379/diff_cpg_fm/data/female_background.RData")
} else{
  female_seq_bg <- makeBackground(motifs = PWMLogn.hg19.MotifDb.Hsap$pwms, bg.seq=femaleDMPseq,
                                              type="logn", algorithm="human",verbose=TRUE)

  save(female_seq_bg, file="/home/og16379/diff_cpg_fm/data/female_background.RData")
}



################################################################################
# comparison
################################################################################

# motifs differentially enriched in the maintained sequence (with lognormal background correction)
if(file.exists("/home/og16379/diff_cpg_fm/data/female_male_enrichmentRep.RData")){
  load("/home/og16379/diff_cpg_fm/data/female_male_enrichmentRep.RData")
} else{
  female_male_enrichment <- motifEnrichment(femaleDMPseq, male_seq_bg)
  save(female_male_enrichment, file="/home/og16379/diff_cpg_fm/data/female_male_enrichmentRep.RData")
}

if(file.exists("/home/og16379/diff_cpg_fm/data/male_female_enrichmentRep.RData")){
  load("/home/og16379/diff_cpg_fm/data/male_female_enrichmentRep.RData")
} else{
  male_female_enrichment <- motifEnrichment(maleDMPseq, female_seq_bg)
  save(male_female_enrichment, file="/home/og16379/diff_cpg_fm/data/male_female_enrichmentRep.RData")
}


# group report
# p value set at 0.05
# get unique motifs and remove UW.Motifs
# write to csv
female_male_enrichment_report <- groupReport(female_male_enrichment)
female_male_enrichment_report_pvalue05 <- female_male_enrichment_report[female_male_enrichment_report$p.value < 0.05]
female_male_enrichment_report_pvalue05 <- female_male_enrichment_report_pvalue05[-grep("UW.Motif.",female_male_enrichment_report_pvalue05$target)]
enrichedInDMPsHMInFemales <- unique(female_male_enrichment_report_pvalue05$target)
female_male_enrichment_report_pvalue05df<-as.data.frame(female_male_enrichment_report_pvalue05)
write.csv(female_male_enrichment_report_pvalue05df,"/home/og16379/diff_cpg_fm/data/femaleEnrichedTFBSRep.csv")

# repeat above for males
male_female_enrichment_report <- groupReport(male_female_enrichment)
male_female_enrichment_report_pvalue05 <- male_female_enrichment_report[male_female_enrichment_report$p.value < 0.05]
male_female_enrichment_report_pvalue05 <- male_female_enrichment_report_pvalue05[-grep("UW.Motif.",male_female_enrichment_report_pvalue05$target)]
enrichedInDMPsHMInMales <- unique(male_female_enrichment_report_pvalue05$target)
male_female_enrichment_report_pvalue05df<-as.data.frame(male_female_enrichment_report_pvalue05)
write.csv(male_female_enrichment_report_pvalue05df,"/home/og16379/diff_cpg_fm/data/maleEnrichedTFBSRep.csv")

# get list of motifs enriched at sites hypermethylated in females only
femaleOnly<-enrichedInDMPsHMInFemales[-which(enrichedInDMPsHMInFemales %in% enrichedInDMPsHMInMales)]
# get list of motifs enriched at sites hypermethylated in males only
maleOnly<-enrichedInDMPsHMInMales[-which(enrichedInDMPsHMInMales %in% enrichedInDMPsHMInFemales)]
# get list of motifs enriched at sites hypermethylated in females and males
overlap<-enrichedInDMPsHMInMales[which(enrichedInDMPsHMInMales %in% enrichedInDMPsHMInFemales)]



############################################
## complete ##
############################################
