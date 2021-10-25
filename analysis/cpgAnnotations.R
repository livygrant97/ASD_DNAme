################################################################
# Characterising Autosomal sex differences between males and females
# Olivia Grant
###############################################################

###############################################################
####  look at genomic annotations ####
###############################################################
#libraries


packages <- c("knitr", "limma", "minfi", "minfi","IlluminaHumanMethylationEPICanno.ilm10b2.hg19","RColorBrewer",
"missMethyl","Gviz","IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
"DMRcate","stringr","bigmelon","BSgenome.Hsapiens.UCSC.hg19","IlluminaHumanMethylationEPICmanifest","org.Hs.eg.db","regioneR","TxDb.Hsapiens.UCSC.hg19.knownGene")

lapply(packages, library, character.only = TRUE)

#load in data
gfile <- openfn.gds('/storage/st05d/Exeter/UnderSocMeth/UScombo/USM_Combo1.gds',readonly=T)

#convert to methylset function
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

# convert gds file to methylumi set file
mSet<-gds2mset(gfile,anno="epic")
mSet<-mapToGenome(mSet,mergeManifest=TRUE)
# get granges
mSet<-granges(mSet)
# load hg19 info
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# get exons
exons<-exons(txdb)
# get promoters
promoters<-promoters(genes)

# load in the rest of the data from UCSC table browser (downloaded to comp)
# transposableElements
transposableElements<-read.table("/home/og16379/diff_cpg_fm/data/repeats.bed")
colnames(transposableElements)<-c("chr","genome","feature","start","end","width" ,"strand","strand2","name","TEtype","x","y","TEtype2","z")
transposableElements<-makeGRangesFromDataFrame(transposableElements,keep.extra.columns=TRUE)
# five prime utrs
five_prime_UTR<-toGRanges("/home/og16379/diff_cpg_fm/data/5UTR.bed")
# three prime utrs
three_prime_UTR<-toGRanges("/home/og16379/diff_cpg_fm/data/3UTR.bed")
# introns
introns<-toGRanges("/home/og16379/diff_cpg_fm/data/UCSC_Introns.tsv")

# enhancer data from enhancer atlas for several cell types in whole blood
enhancerHSC<-toGRanges("/home/og16379/diff_cpg_fm/data/HSC.bed")
enhancerMonocyteCD14<-toGRanges("/home/og16379/diff_cpg_fm/data/CD14+_monocyte.bed")
enhancerMonocyte<-toGRanges("/home/og16379/diff_cpg_fm/data/Monocyte.bed")

#######################
# get overlaps
#######################

five_prime_UTRov<-findOverlaps(five_prime_UTR,mSet)
five_prime_UTRovx<-unique(subjectHits(five_prime_UTRov))
cgin5prime<-DataFrame(cg=names(mSet)[five_prime_UTRovx])
cgin5prime$feature="5' UTR"

three_prime_UTRov<-findOverlaps(three_prime_UTR,mSet)
three_prime_UTRovx<-unique(subjectHits(three_prime_UTRov))
cgin3prime<-DataFrame(cg=names(mSet)[three_prime_UTRovx])
cgin3prime$feature="3' UTR"

promoterov<-findOverlaps(promoters,mSet)
promoterovx<-unique(subjectHits(promoterov))
cginpromoter<-DataFrame(cg=names(mSet)[promoterovx])
cginpromoter$feature="Promoter"

transposableElementsov<-findOverlaps(transposableElements,mSet)
transposableElementsovx<-unique(subjectHits(transposableElementsov))
cginte<-DataFrame(cg=names(mSet)[transposableElementsovx])
cginte$feature="TE"

transposableElementsov<-findOverlaps(transposableElements,mSet)
transposableElementsovx<-unique(subjectHits(transposableElementsov))
cginte<-DataFrame(cg=names(mSet)[transposableElementsovx])
cginte$feature="TE"

intronsov<-findOverlaps(introns,mSet)
intronsovx<-unique(subjectHits(intronsov))
cginintrons<-DataFrame(cg=names(mSet)[intronsovx])
cginintrons$feature="Introns"

genesov<-findOverlaps(genes,mSet)
genesovx<-unique(subjectHits(genesov))
cgingenes<-DataFrame(cg=names(mSet)[genesovx])
cgingenes$feature="Genes"

cpgislandsov<-findOverlaps(cpgIslands,mSet)
cpgislandsovx<-unique(subjectHits(cpgislandsov))
cgincpgislands<-DataFrame(cg=names(mSet)[cpgislandsovx])
cgincpgislands$feature="CpG Islands"

exonsov<-findOverlaps(exons,mSet)
exonsovx<-unique(subjectHits(exonsov))
cginexons<-DataFrame(cg=names(mSet)[exonsovx])
cginexons$feature="Exons"

enhancerHSCov<-findOverlaps(enhancerHSC,mSet)
enhancerHSCovx<-unique(subjectHits(enhancerHSCov))
cginenhancerHSC<-DataFrame(cg=names(mSet)[enhancerHSCovx])
cginenhancerHSC$feature="enhancer"

enhancerMonocyteCD14ov<-findOverlaps(enhancerMonocyteCD14,mSet)
enhancerMonocyteCD14ovx<-unique(subjectHits(enhancerMonocyteCD14ov))
cginenhancerCD14<-DataFrame(cg=names(mSet)[enhancerMonocyteCD14ovx])
cginenhancerCD14$feature="enhancer"

enhancerMonocyteov<-findOverlaps(enhancerMonocyte,mSet)
enhancerMonocyteovx<-unique(subjectHits(enhancerMonocyteov))
cginenhancerMon<-DataFrame(cg=names(mSet)[enhancerMonocyteovx])
cginenhancerMon$feature="enhancer"

cgAnnotations<-rbind(cginpromoter,cginte,cginintrons,cgin3prime,cgin5prime,cginexons,cginenhancerHSC,cginenhancerCD14,cginenhancerMon)
