

  ################################################################
  # Characterising Autosomal sex differences between males and females
  # Olivia Grant
  ###############################################################

  # below is all the code to produce the figures shown in figure 3A-3D
  # previous code must have been run in order to produce these figures


  ################################################################################
  # load libraries
  ################################################################################

  if(!require("bigmelon", character.only = TRUE)){
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    BiocManager::install("bigmelon")
  }
  library(bigmelon)

  if(!require("Gviz", character.only = TRUE)){
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    BiocManager::install("Gviz")
  }
  library(Gviz)

  if(!require("IlluminaHumanMethylationEPICanno.ilm10b2.hg19", character.only = TRUE)){
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
  }
  library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)


  ################################################################################
  # prepare data for plotting
  ################################################################################

  # splitting by sex
  # convert to df
  dasen_autosome<-dasen_autosome[,!(colnames(dasen_autosome) %in% outlier)]
  dasendf<-as.data.frame(dasen_autosome)
  # extract female samples
  femaleBetas1<-dasendf[,info$nsex=="F"]
  # extract male samples
  maleBetas1<-dasendf[,info$nsex=="M"]

  femaleBetas1<-ecc[,info$nsex=="F"]
  # extract male samples
  maleBetas1<-ecc[,info$nsex=="M"]

  # now we will make a granges object for the epic data set
  # re code the gds2mset function from big melon with right annotation

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


  # convert gds to mset object
  gfile1<-gds2mset(gfile)

  # map it to the genome
  msetMapped<-mapToGenome(gfile1)

  # convert that to a g range object
  msetMapped<-granges(msetMapped)

  # subset to get saDMPs
  femaleGRanges<-subset(msetMapped,msetMapped@ranges@NAMES %in% repFemalesSigCpG2$Row.names)
  maleGRanges<-subset(msetMapped,msetMapped@ranges@NAMES %in% repMalesSigCpG2$Row.names)


  #DMRsexSig<-get(load("/home/og16379/diff_cpg_fm/data/DMRsex.Rdata"))

  # load in extra libraries
  library(rtracklayer)
  library(RColorBrewer)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)

  #get txdb object
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  #set genome
  gen <- "hg19"
  # the index of the DMR that we will plot (simply change this to dmr of interest for each plot, here I show how to plot dmr 1)
  dmrIndex <- 4
  #set col
  pal <- brewer.pal(2,"Dark2")
  #set info for plotting
  chrom <- as.character(seqnames(DMRsexSig[dmrIndex]))
  start <- as.numeric(start(DMRsexSig[dmrIndex]))
  end <- as.numeric(end(DMRsexSig[dmrIndex]))
  #add extra 25% to plot, change as you see fit for each plot
  minbase <- start - (0.25*(end-start))
  maxbase <- end + (0.25*(end-start))

  # genome axis track
  gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
  #gene region track
  grTrack <- GeneRegionTrack(txdb, genome = gen,
                            chromosome = chrom, name = "Genes",
                            start=start, to=end,collapseTranscripts="meta",col.title="black",
                            showID=FALSE,stacking="squish")

  # i set ideoTrack manually using cytoband data from UCSC
  data<-read.table("/home/og16379/diff_cpg_fm/data/cytoBand.txt",header=F,sep="\t")
  colnames(data) <-c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain')
  ideoTrack <- IdeogramTrack(genome="hg19",chromosome = chrom , bands=data, from=start,to=end)

  data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
  annotationTable <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  annotationTable<-as.data.frame(annotationTable)
  annotationTableOrd <- annotationTable[order(annotationTable$chr,annotationTable$pos),]
  head(annotationTableOrd)

  # make sure all in same order for plotting
  femaleBetasOrd <- femaleBetas1[match(annotationTableOrd$Name,rownames(femaleBetas1)),]
  head(femaleBetasOrd)
  maleBetasOrd <- maleBetas1[match(annotationTableOrd$Name,rownames(maleBetas1)),]
  head(maleBetasOrd)

  #dmr track

  dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR",
  chromosome=chrom,fill=cbbPalette[4],col.title="black")

  # we want to plot with m values instead of beta values so convert beta to m value
  femaleM<-beta2m(femaleBetasOrd)
  maleM<-beta2m(maleBetasOrd)

  # set up female and male track data
  femalecpgData <- GRanges(seqnames=Rle(annotationTableOrd$chr),
  ranges=IRanges(start=annotationTableOrd$pos, end=annotationTableOrd$pos),
  strand=Rle(rep("*",nrow(annotationTableOrd))),
  betas=femaleM)

  malecpgData <- GRanges(seqnames=Rle(annotationTableOrd$chr),
  ranges=IRanges(start=annotationTableOrd$pos, end=annotationTableOrd$pos),
  strand=Rle(rep("*",nrow(annotationTableOrd))),
  betas=maleM)

  femalecpgData <- subsetByOverlaps(femalecpgData, DMRsexSig[dmrIndex])
  malecpgData <- subsetByOverlaps(malecpgData, DMRsexSig[dmrIndex])

  maleTrack <- DataTrack(range=malecpgData,genome = gen,  col=maleColour,
     chromsome=chrom,
     type=c("a","p"),  name="Methylation values",
     background.panel="white",col.title="black",cex=0.1,showAxis=TRUE)

  femaleTrack <- DataTrack(range=femalecpgData,genome = gen,  col=femaleColour,
     chromsome=chrom,
     type=c("a","p"),  name="Methylation values",
     background.panel="white",col.title="black",cex=0.1,showAxis=TRUE)

  #overlay tracks for plotting
  ot<-OverlayTrack(trackList=list(femaleTrack,maleTrack),name="male vs Female")


  #set sizes of tracks
  sizes <- c(1,1,6,1,1) # set up the relative sizes of the tracks

  #plot dmr
  pdf("/home/og16379/diff_cpg_fm/pdf/dmrTrack4Rep.pdf")
  plotTracks(list(ideoTrack,gTrack,ot,dmrTrack,grTrack),from=minbase,to=maxbase,showTitle=TRUE,sizes=sizes)
  dev.off()

  # repeat for each DMR of interest

  #########################
