################################################################
# Characterising Autosomal sex differences between males and females
# Olivia Grant
###############################################################

################################################################
# calculating DMRs for discovery data set
###############################################################


# outliers
outlier <- c('200611820013_R08C01', '200603220075_R08C01', '200864580018_R02C01', '200611820020_R07C01')
# calculate cell type prop estimates
cellCounts<-estimateCellCounts.gds(gfile)
# remove outliers
cellCounts<-cellCounts[!(rownames(cellCounts) %in% outlier),]
# create pheno data
info <- pData(gfile)[, c('barcode', 'nsex','confage')]
info$nsex <- gsub('2', 'F', info$nsex)
info$nsex <- gsub('1', 'M', info$nsex)
info <- info[!(info$barcode %in% outlier), ]
dasen_autosome <- dasen_autosome[, info$barcode]
phenoGroup<-cbind(info,cellCounts)
rownames(phenoGroup)<-colnames(dasen_autosome)
designSex<-model.matrix(~nsex+confage+CD8T+CD4T+NK+Bcell+Mono+Gran,data=phenoGroup)

dmrSex<-cpg.annotate(object=dasen_autosome,datatype="array",what="B",analysis.type="differential",design=designSex ,coef="nsexM")
dmrSex<-dmrcate(dmrSex)
DMRsexSig<-extractRanges(dmrSex)
save(DMRsexSig,file="/home/og16379/diff_cpg_fm/data/DMRsexFirstRound.Rdata")
save(dmrSex,file="/home/og16379/diff_cpg_fm/data/DMRsexcateFirstRound.Rdata")

DMRdf<-as.data.frame(DMRsexSig)
DMRsig<-subset(DMRdf,  (minfdr < 0.05 & maxbetafc > 0.05 | maxbetafc < -0.05) & no.cpgs > 5)
DMRsig<- DMRsig[!DMRsig$seqnames %in% c("chrX","chrY"),]
DMRsexSig<-makeGRangesFromDataFrame(DMRsig,keep.extra.columns=TRUE)
save(DMRsexSig,file="/home/og16379/diff_cpg_fm/data/DMRsexSignificantFirstRound.Rdata")
write.csv(DMRsexSig, file="/home/og16379/diff_cpg_fm/data/saDMRAdditionalFileFirstRound.csv")


################################################################
# calculating DMRs for validation data set
###############################################################
load("//home/og16379/diff_cpg_fm/data/dasen_autosome2round.Rdata")

# outliers
outlier2 <- c('203991410050_R06C01', '203991410050_R07C01', '204026590012_R07C01', '203998240051_R07C01',
    '204026590078_R04C01', '203994670097_R08C01', '203960330159_R08C01', '204022160070_R03C01',
    '203994670097_R02C01', '203991470090_R05C01', '204026590012_R06C01')

# calculate cell type prop estimates
cellCounts2<-estimateCellCounts.gds(gfile_2)
# remove outliers
cellCounts2<-cellCounts2[!(rownames(cellCounts2) %in% outlier2),]
# create pheno data
info2 <- pData(gfile_2)[, c('barcode','Sex','Age')]
info2$Sex <- gsub('2', 'F', info2$Sex)
info2$Sex <- gsub('1', 'M', info2$Sex)
info2 <- info2[!(info2$barcode %in% outlier2), ]
dasen_autosome2 <- dasen_autosome2[, info2$barcode]
phenoGroup2<-cbind(info2,cellCounts2)
rownames(phenoGroup2)<-colnames(dasen_autosome2)
designSex2<-model.matrix(~Sex+Age+CD8T+CD4T+NK+Bcell+Mono+Gran,data=phenoGroup2)

dmrSex2<-cpg.annotate(object=dasen_autosome2,datatype="array",what="B",analysis.type="differential",design=designSex2 ,coef="SexM")
dmrSex2<-dmrcate(dmrSex2)
DMRsexSig2<-extractRanges(dmrSex2)
save(DMRsexSig2,file="/home/og16379/diff_cpg_fm/data/DMRsexSecondRound.Rdata")
save(dmrSex2,file="/home/og16379/diff_cpg_fm/data/DMRsexcateSecondRound.Rdata")

DMRdf2<-as.data.frame(DMRsexSig2)
DMRsig2<-subset(DMRdf2,  (minfdr < 0.05 & maxbetafc > 0.05 | maxbetafc < -0.05) & no.cpgs > 5)
DMRsig2<- DMRsig2[!DMRsig2$seqnames %in% c("chrX","chrY"),]
DMRsexSig2<-makeGRangesFromDataFrame(DMRsig2,keep.extra.columns=TRUE)
save(DMRsexSig2,file="/home/og16379/diff_cpg_fm/data/DMRsexSignificantSecondRound.Rdata")
write.csv(DMRsexSig2, file="/home/og16379/diff_cpg_fm/data/saDMRAdditionalFileSecondRound.csv")


################################################################
# see which DMRs are replicated across data sets
###############################################################

# subset by overlaps and order by FDR
significantSaDMRs<- subsetByOverlaps(DMRsexSig,DMRsexSig2)
sDf <- significantSaDMRs[order(significantSaDMRs$minfdr,decreasing=FALSE),]
write.csv(sDf, file="/home/og16379/diff_cpg_fm/data/saDMRAdditionalFile.csv")


################################################################
# pathway analysis on saDMRs
###############################################################

anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

  cpgs <- GenomicRanges::GRanges(seqnames = anno$chr,
                  ranges = IRanges::IRanges(start = anno$pos,
                                   end = anno$pos),
                  strand = anno$strand,
                  name = anno$Name)

  overlaps <- GenomicRanges::findOverlaps(cpgs,sDf)
  sig.cpg <- cpgs$name[from(overlaps)]

  resultKEGG <- gometh(sig.cpg=sig.cpg, all.cpg=all.cpg, collection="KEGG")
  # nothing significant
  resultGO <- gometh(sig.cpg=sig.cpg, all.cpg=all.cpg, collection="GO")
 # nothing significant


