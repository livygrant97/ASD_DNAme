################################################################
# Characterising Autosomal sex differences between males and females
# Olivia Grant
###############################################################

#below is all the code to produce the figures shown in figure 1A-1E
#previous code must have been run in order to produce these figures

###############################################################
####  Generate Figure 1A  ####
###############################################################


# subset annotation table
annotationTable1<-subset(annotationTable,rownames(annotationTable) %in% sigCpG_Autosomes_Cleaned$Row.names)
sigCpGs<-subset(sigCpG_Autosomes_Cleaned,Row.names %in% rownames(annotationTable1))

#create a data frame for plotting
asd_Manhattan<- data.frame('CHR' = as.character(annotationTable1$chr),
                   'BP' = as.numeric(annotationTable1$pos),
                   'P' = as.numeric(abs(sigCpGs$t)),
                   "GENE"= as.character(sigCpGs$gene),
                   "SNP"=as.character(sigCpGs$Row.names),
                     stringsAsFactors = F
                   )

# Replace names of X and Y to 23 and 24
seqlevelsStyle(asd_Manhattan $CHR)<-"NCBI"
asd_Manhattan$CHR <- as.numeric(asd_Manhattan $CHR)
#we are highlighting diff chromosomes and sadmps on the chromosomes
is.even <- function(x) x %% 2 == 0
is.odd <- function(x) x %% 2 != 0
hlight<-asd_Manhattan[asd_Manhattan$SNP %in% replicatedResults$Row.names,]
hlightgroup1<- hlight[is.even(hlight$CHR),]
hlightgroup2<- hlight[is.odd(hlight$CHR),]
grp1<-hlightgroup1$SNP
grp2<-hlightgroup2$SNP


# adapt manhattan function from qqman
manhattanLivsFunction<-function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10",
"gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05),annotatePval=NULL, annotateTop=TRUE,
genomewideline = -log10(5e-08), highlight1 = NULL, highlight2 = NULL, logp = TRUE,
...)
{
CHR = BP = P = index = NULL
if (!(chr %in% names(x)))
stop(paste("Column", chr, "not found!"))
if (!(bp %in% names(x)))
stop(paste("Column", bp, "not found!"))
if (!(p %in% names(x)))
stop(paste("Column", p, "not found!"))
if (!(snp %in% names(x)))
warning(paste("No SNP column found. OK unless you're trying to highlight."))
if (!is.numeric(x[[chr]]))
stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
if (!is.numeric(x[[bp]]))
stop(paste(bp, "column should be numeric."))
if (!is.numeric(x[[p]]))
stop(paste(p, "column should be numeric."))
d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
if (!is.null(x[[snp]]))
d = transform(d, SNP = x[[snp]])
d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
d <- d[order(d$CHR, d$BP), ]
if (logp) {
d$logp <- -log10(d$P)
}
else {
d$logp <- d$P
}
d$pos = NA
d$index = NA
ind = 0
for (i in unique(d$CHR)) {
ind = ind + 1
d[d$CHR == i, ]$index = ind
}
nchr = length(unique(d$CHR))
if (nchr == 1) {
options(scipen = 999)
d$pos = d$BP/1e+06
ticks = floor(length(d$pos))/2 + 1
xlabel = paste("Chromosome", unique(d$CHR), "position(Mb)")
labs = ticks
}
else {
lastbase = 0
ticks = NULL
for (i in unique(d$index)) {
if (i == 1) {
d[d$index == i, ]$pos = d[d$index == i, ]$BP
}
else {
lastbase = lastbase + tail(subset(d, index ==
i - 1)$BP, 1)
d[d$index == i, ]$pos = d[d$index == i, ]$BP +
lastbase
}
ticks = c(ticks, (min(d[d$CHR == i, ]$pos) + max(d[d$CHR ==
i, ]$pos))/2 + 1)
}
xlabel = "Chromosome"
labs <- unique(d$CHR)
}
xmax = ceiling(max(d$pos) * 1.03)
xmin = floor(max(d$pos) * -0.03)
def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i",
las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0,140), xlab = xlabel, ylab = expression(T-Value))
dotargs <- list(...)
do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in%
names(dotargs)]))
if (!is.null(chrlabs)) {
if (is.character(chrlabs)) {
if (length(chrlabs) == length(labs)) {
labs <- chrlabs
}
else {
warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
}
}
else {
warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
}
}
if (nchr == 1) {
axis(1, ...)
}
else {
axis(1, at = ticks, labels = labs, ...)
}
col = rep(col, max(d$CHR))
if (nchr == 1) {
with(d, points(pos, logp, pch = 20, col = col[1], ...))
}
else {
icol = 1
for (i in unique(d$index)) {
with(d[d$index == unique(d$index)[i], ], points(pos,
logp, col = col[icol], pch = 20, ...))
icol = icol + 1
}
}
if (suggestiveline)
abline(h = suggestiveline, col = "blue")
if (genomewideline)
abline(h = genomewideline, col = "grey")
if (!is.null(highlight1)) {
if (any(!(highlight1 %in% d$SNP)))
warning("You're trying to highlight1 SNPs that don't exist in your results.")
d.highlight1 = d[which(d$SNP %in% highlight1), ]
with(d.highlight1, points(pos, logp, col = "#0072B2", pch = 20,
...))
}
if (!is.null(highlight2)) {
if (any(!(highlight2 %in% d$SNP)))
warning("You're trying to highlight2 SNPs that don't exist in your results.")
d.highlight2 = d[which(d$SNP %in% highlight2), ]
with(d.highlight2, points(pos, logp, col = "#D55E00", pch = 20,
...))
}

if (!is.null(annotatePval)) {
    # extract top SNPs at given p-val
    if (logp) {
        topHits = subset(d, P <= annotatePval)
    } else
        topHits = subset(d, P >= annotatePval)
    par(xpd = TRUE)
    # annotate these SNPs
    if (annotateTop == FALSE) {
      if (logp) {
          with(subset(d, P <= annotatePval),
               textxy(pos, -log10(P), offset=0.625, labs = topHits$SNP, cex = 0.10), ...)
      } else
          with(subset(d, P >= annotatePval),
               textxy(pos, P, offset = 0.625, labs = topHits$SNP, cex = 0.10), ...)
    }
    else {
        # could try alternative, annotate top SNP of each sig chr
        topHits <- topHits[order(topHits$P),]
        topSNPs <- NULL

        for (i in unique(topHits$CHR)) {

            chrSNPs <- topHits[topHits$CHR == i,]
            topSNPs <- rbind(topSNPs, chrSNPs[1,])
                }
        if (logp ){
            textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, labs = topSNPs$SNP, cex = 0.5, ...)
        } else
          textxy(topSNPs$pos, topSNPs$P, offset = 0.4, labs = topSNPs$SNP, cex = 0.5, ...)
    }
}
par(xpd = FALSE)
}


# plotting using manhattan
dmpManhattan <- function(asd_Manhattan){
library(png)
library(calibrate)

           fold.range <- c(0,449212500)
            prob.range <- c(0,140)
            plot(
                    x = 0, y = 0, xlim = fold.range, ylim = prob.range,
                    xlab = "", ylab = "",
                    xaxt="n", type = "n", cex.axis = 1.0, cex.lab = 1.0,yaxs="i",xaxs="i")

            coords <- par("usr")
            gx <- grconvertX(coords[1:2], "user", "in")
            gy <- grconvertY(coords[3:4], "user", "in")
            width <- max(gx) - min(gx)
            height <- max(gy) - min(gy)

            tmp <- tempfile()
            png(tmp, width = width, height = height, units = "in", res = 500, bg = "transparent")

            plot.new()
            par(mar = c(0.0, 0.0, 0.0, 0.0))

           manhattanLivsFunction(asd_Manhattan,highlight=grp1,highlight2=grp2,logp=FALSE,
           col=c(cbPalette[2],cbPalette[3]),suggestiveline = F, genomewideline = 5,annotatePval=40)
           dev.off()
            panel <- readPNG(tmp)
            rasterImage(panel, coords[1], coords[3], coords[2], coords[4])

            mtext(text = "abs(T-value)", side = 2, line = 2.5, cex = 1)
            mtext(text="Chromosome",side=1,line=1.5,cex=1)

   }


    png("/home/og16379/diff_cpg_fm/png/Manhattan.png",width=3000,height=3000,units="px",res=400)
    manhattanLivsFunction(asd_Manhattan,highlight=grp1,highlight2=grp2,logp=FALSE,
    col=c(cbPalette[2],cbPalette[3]),suggestiveline = F, genomewideline = 5,annotatePval=40)
    dev.off()

###############################################################
####  Generate Figure 1B  ####
###############################################################
# note: we need to filter differently here - we don't remove probes that
# dont pass the FDR test because we want to plot the NS sites too
# detirmine what sex methylation is higher in
# add a column called higher met in which indicates which sex met is higher in

dmp$Higher_Met_In<-"Not significant"
dmp$Higher_Met_In [dmp$Row.names %in% repFemalesSigCpG2$Row.names] <- "Female"
dmp$Higher_Met_In [dmp$Row.names %in% repMalesSigCpG2$Row.names] <- "Male"
# same as below
data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
annotationTable <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
annotationTable<-as.data.frame(annotationTable)
#get probes on x and y chromosome: 19681 probes
sexProbes<-as.character(annotationTable$Name[annotationTable$chr %in% c("chrX","chrY")])
#remove probes not on autosomes
all_Autosomes<- dmp[!dmp$Row.names %in% sexProbes,]
# read in annotations for SNP's and cross hybridizing probes and remove them from analysis
# https://zwdzwd.github.io/InfiniumAnnotation
SNPAnnotations<-read.delim("/home/og16379/diff_cpg_fm/data/EPIC.hg19.manifest.tsv")
SNPAssociatedCpG<- SNPAnnotations %>% dplyr::filter(MASK_general==TRUE)
all_Autosomes_Cleaned<- all_Autosomes[!all_Autosomes$Row.names %in%  SNPAssociatedCpG$probeID,]
# downloaded the list of cross hybridizing probes from here https://pubmed.ncbi.nlm.nih.gov/23314698/
chProbes <- read.table("/home/og16379/diff_cpg_fm/data/chProbes.txt")
all_Autosomes_Cleaned<-all_Autosomes_Cleaned[!(all_Autosomes_Cleaned$Row.names
  %in% chProbes$V1),]

#get the top hits so we can add labels to the plot
topHits<-rbind(repFemalesSigCpG2[1:8,],repMalesSigCpG2[1:8,])

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
maleColour<-cbbPalette[7]
femaleColour<-cbbPalette[6]

all_Autosomes_Cleaned$Colour<-"#999999"
all_Autosomes_Cleaned$Colour[all_Autosomes_Cleaned$Row.names %in% repFemalesSigCpG2$Row.names] <- femaleColour
all_Autosomes_Cleaned$Colour[all_Autosomes_Cleaned$Row.names %in% repMalesSigCpG2$Row.names] <- maleColour


dmpVolcano <- function(all_Autosomes_Cleaned,deltaBetaLimit,tLimit){

            library(maptools)
            fold.range <- c(-0.32,0.32)
            prob.range <- range(log10(abs(all_Autosomes_Cleaned[,'t'])))

            plot(
                    x = 0, y = 0, xlim = fold.range, ylim = prob.range,
                    xlab = "", ylab = "",
                    axes = TRUE, type = "n", cex.axis = 1.0, cex.lab = 1.0, cex.main = 2.0)

            coords <- par("usr")
            gx <- grconvertX(coords[1:2], "user", "in")
            gy <- grconvertY(coords[3:4], "user", "in")
            width <- max(gx) - min(gx)
            height <- max(gy) - min(gy)

            tmp <- tempfile()
            png(tmp, width = width, height = height, units = "in", res = 500, bg = "transparent")

                    plot.new()
                    par(mar = c(0.0, 0.0, 0.0, 0.0))
                    plot(
                            x = 0, y = 0, ,xlim = fold.range, ylim = prob.range,
                            xlab = "", ylab = "", #main = title,
                            axes = FALSE, type = "n")


                    x.data = all_Autosomes_Cleaned[,'deltaBeta']
                    y.data = log10(abs(all_Autosomes_Cleaned[,'t']))
                    y.data[y.data == Inf] = 0
                    y.data[is.na(y.data)] = 0
                    points(x = x.data, y = y.data,col=all_Autosomes_Cleaned$Colour,cex = 0.8,pch=16)

                   abline(v = -deltaBetaLimit,col="grey")
                    abline(v = deltaBetaLimit, col="grey")
                    abline(h = log10(abs(tLimit)), col="grey")
                    pointLabel(topHits[,'deltaBeta'],log10(abs(topHits[,'t'])),labels=paste("",topHits$gene,"",sep=""),cex=0.5)

            dev.off()
            panel <- readPNG(tmp)
            rasterImage(panel, coords[1], coords[3], coords[2], coords[4])
            mtext(text = expression(Delta*"Beta male-female"), side = 1, line = 4, cex = 1)
            mtext(text = "T-value", side = 2, line = 2.5, cex = 1)

    }

pdf("/home/og16379/diff_cpg_fm/pdf/volcanoPlotRep.pdf")
    dmpVolcano(all_Autosomes_Cleaned,0.05,5)
dev.off()



###############################################################
####  Generate Figure 1C ####
###############################################################

# subset the beta matrix for only the significant saDMPs
library(factoextra)
betas<- read.gdsn(index.gdsn(gfile,"betas"))
rownames(betas)<-read.gdsn(index.gdsn(gfile,"fData/Name"))
colnames(betas)<-read.gdsn(index.gdsn(gfile,"pData/barcode"))

info <- info[!(info$barcode %in% outlier), ]
betas <- betas[, info$barcode]

sigProbeBetas <- subset(betas,rownames(betas) %in% replicatedResults$Row.names)
# perform principal components analysis
pca <- prcomp(t(sigProbeBetas)) ##### samples in row, probes in column !!!
#### print the explained variance in the first five PCs
summary(pca)$importance[, 1:5]
#set f and m to female and male
grp<-as.factor(info$nsex)
grp <- gsub('F', 'Females', grp)
grp<- gsub('M', 'Males', grp)
#plot figure 1b
pdf("/home/og16379/diff_cpg_fm/pdf/pcaPlotRep.pdf",width=6,height=5)
fviz_pca_ind(pca, geom.ind="point",pointshape=19,habillage=grp,palette=c("#0072B2", "#D55E00"),legend.title="Sex",addEllipses=FALSE,title="",alpha.ind=0.5)
dev.off()

###############################################################
####  Generate Figure 1D ####
###############################################################

# now we will calculate how many saDMPs per gene to see distribution of saDMPs in diff genes
geneDuplicated<- as.data.frame(replicatedResults %>% group_by(gene) %>% tally())
geneDuplicated<-geneDuplicated[-1,]
pdf("/home/og16379/diff_cpg_fm/pdf/saDMPperGenePlotRep.pdf",width=6)
barplot(table(geneDuplicated$n),xlab="Number of saDMPs per gene",ylab="Frequency",col=cbbPalette[4],xlim=c(0,9))
dev.off()

###############################################################
####  Generate Figure 1E ####
###############################################################

#libraries

packages <- c("knitr", "limma", "minfi", "minfi","IlluminaHumanMethylationEPICanno.ilm10b2.hg19","RColorBrewer",
"missMethyl","Gviz","IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
"DMRcate","stringr","bigmelon","BSgenome.Hsapiens.UCSC.hg19","IlluminaHumanMethylationEPICmanifest","org.Hs.eg.db","regioneR","TxDb.Hsapiens.UCSC.hg19.knownGene")

lapply(packages, library, character.only = TRUE)


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
# get genes
genes<-genes(txdb)
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

df<-data.frame("CG"=as.character(rownames(annotationTable)),"Feature"="other")
df$CG<-as.character(df$CG)
df$Feature<-as.character(df$Feature)
df <- df[match(mSet@ranges@NAMES,df$CG),]

df$Feature[df$CG %in% cginexons$cg] <- "Exon"
df$Feature[df$CG %in% cginintrons$cg] <- "Introns"
df$Feature[df$CG %in% cgin5prime$cg] <- "5' UTR"
df$Feature[df$CG %in% cgin3prime$cg] <- "3' UTR"
df$Feature[df$CG %in% cginenhancerHSC$cg] <- "Enhancers"
df$Feature[df$CG %in% cginenhancerCD14$cg] <- "Enhancers"
df$Feature[df$CG %in% cginenhancerMon$cg] <- "Enhancers"
df$Feature[df$CG %in% cginpromoter$cg] <- "Promoters"
df$Feature[df$CG %in% cginte$cg] <- "TE's"

df$Group<- "Background"
df$Group[df$CG %in% repMalesSigCpG$Row.names] <- "males"
df$Group[df$CG %in% repFemalesSigCpG$Row.names] <- "females"

results_matrix <- matrix(0, 4, 8)
colnames(results_matrix) <- c("Other", "Exons", "Introns", "5' UTR", "3' UTR", "Enhancers", "Promoters","TE")
rownames(results_matrix) <- c("All Autosomes", "all saDMPs", "Hypermethylated in females", "Hypermethylated in males")



results_matrix[1,1] <- length(which(df$Group == "Background" & df$Feature== "other"))
results_matrix[1,2] <- length(which(df$Group == "Background" & df$Feature == "Exon"))
results_matrix[1,3] <- length(which(df$Group == "Background" & df$Feature == "Introns"))
results_matrix[1,4] <- length(which(df$Group == "Background" & df$Feature == "5' UTR"))
results_matrix[1,5] <- length(which(df$Group == "Background" & df$Feature == "3' UTR"))
results_matrix[1,6] <- length(which(df$Group == "Background" & df$Feature == "Enhancers"))
results_matrix[1,7] <- length(which(df$Group == "Background" & df$Feature == "Promoters"))
results_matrix[1,8] <- length(which(df$Group == "Background" & df$Feature== "TE's"))

results_matrix[2,1] <- length(which(df$Group == "females" | df$Group == "males" & df$Feature== "other"))
results_matrix[2,2] <- length(which(df$Group == "females" | df$Group == "males" & df$Feature == "Exon"))
results_matrix[2,3] <- length(which(df$Group == "females" | df$Group == "males" & df$Feature == "Introns"))
results_matrix[2,4] <- length(which(df$Group == "females" | df$Group == "males" & df$Feature == "5' UTR"))
results_matrix[2,5] <- length(which(df$Group == "females" | df$Group == "males" &  df$Feature == "3' UTR"))
results_matrix[2,6] <- length(which(df$Group == "females" | df$Group == "males" & df$Feature == "Enhancers"))
results_matrix[2,7] <- length(which(df$Group == "females" | df$Group == "males" & df$Feature == "Promoters"))
results_matrix[2,8] <- length(which(df$Group == "females" | df$Group == "males" & df$Feature== "TE's"))

results_matrix[3,1] <- length(which(df$Group == "females" & df$Feature== "other"))
results_matrix[3,2] <- length(which(df$Group == "females" & df$Feature == "Exon"))
results_matrix[3,3] <- length(which(df$Group == "females" & df$Feature == "Introns"))
results_matrix[3,4] <- length(which(df$Group == "females" & df$Feature == "5' UTR"))
results_matrix[3,5] <- length(which(df$Group == "females" & df$Feature == "3' UTR"))
results_matrix[3,6] <- length(which(df$Group == "females" & df$Feature == "Enhancers"))
results_matrix[3,7] <- length(which(df$Group == "females" & df$Feature == "Promoters"))
results_matrix[3,8] <- length(which(df$Group == "females" & df$Feature== "TE's"))

results_matrix[4,1] <- length(which(df$Group == "males" & df$Feature== "other"))
results_matrix[4,2] <- length(which(df$Group == "males" & df$Feature == "Exon"))
results_matrix[4,3] <- length(which(df$Group == "males" & df$Feature == "Introns"))
results_matrix[4,4] <- length(which(df$Group == "males" & df$Feature == "5' UTR"))
results_matrix[4,5] <- length(which(df$Group == "males" & df$Feature == "3' UTR"))
results_matrix[4,6] <- length(which(df$Group == "males" & df$Feature == "Enhancers"))
results_matrix[4,7] <- length(which(df$Group == "males" & df$Feature == "Promoters"))
results_matrix[4,8] <- length(which(df$Group == "males" & df$Feature== "TE's"))



hmap <- results_matrix
results_matrix <- t(results_matrix)
results_matrix_aves <- results_matrix

results_matrix_aves[,1] <- results_matrix[,1]/sum(results_matrix[,1])
results_matrix_aves[,2] <- results_matrix[,2]/sum(results_matrix[,2])
results_matrix_aves[,3] <- results_matrix[,3]/sum(results_matrix[,3])
results_matrix_aves[,4] <- results_matrix[,4]/sum(results_matrix[,4])

pdf("/home/og16379/diff_cpg_fm/pdf/stackedBarplotRep.pdf")
par(mar=c(4,4,4,10.5))
par(xpd = NA)
barplot(results_matrix_aves, col=cbPalette, yaxt = "none",cex.names=.4)
axis(2, seq(0, 1, 0.1), las=2,labels = paste0(seq(0, 100, 10), "%"))
legend(5,0.75,rownames(results_matrix_aves),fill=cbPalette, bty = "n")
dev.off()



hmap_aves <- hmap


hmap_aves[1,] <- hmap[1,]/sum(hmap[1,])
hmap_aves[2,] <- hmap[2,]/sum(hmap[2,])
hmap_aves[3,] <- hmap[3,]/sum(hmap[3,])
hmap_aves[4,] <- hmap[4,]/sum(hmap[4,])


hmap_comp <- hmap_aves[2:4,]
hmap_comp[1,] <- log2(hmap_comp[1,]/hmap_aves[1,])
hmap_comp[2,] <- log2(hmap_comp[2,]/hmap_aves[1,])
hmap_comp[3,] <- log2(hmap_comp[3,]/hmap_aves[1,])

hmap_comp[hmap_comp > 2] <- 2
hmap_comp[hmap_comp < -2] <- -2


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
custom_at <- seq(-2.5, 2.5, by = (5/30))
cols_contrast <- rev(c(colorRampPalette(c(cbbPalette[7], "white", cbbPalette[6]))(30)))


rownames(hmap_comp)<-c("saDMPs","HM in females","HM in males")
library(lattice)

pdf("/home/og16379/diff_cpg_fm/pdf/heatmap.pdf")
  levelplot(t(hmap_comp),
    at = custom_at,
    xlab = "log2 (obs/exp)",ylab="",
    scales=list(x=list(rot=90)),
    col.regions = cols_contrast)
dev.off()

###############################################################
####  Generate Figure 1F ####
###############################################################


df2<-as.data.frame(table(annotationTable$Relation_to_Island))
colnames(df2)<-c("Relation","Freq")
df2$Group<-"Background"

femaleFeature<-subset(annotationTable,Name %in% repFemalesSigCpG$Row.names)
femaleFeature<-as.data.frame(table(femaleFeature$Relation_to_Island))

maleFeature<-subset(annotationTable,Name %in% repMalesSigCpG$Row.names)
maleFeature<-as.data.frame(table(maleFeature$Relation_to_Island))

saDMPFeature<-subset(annotationTable,Name %in% replicatedResults$Row.names)
saDMPFeature<-as.data.frame(table(saDMPFeature$Relation_to_Island))

backgroundFeature<-as.data.frame(table(annotationTable$Relation_to_Island))


results_matrix2 <- matrix(0, 4, 4)
colnames(results_matrix2) <- c("Island","Shelf","Open sea","Shore")
rownames(results_matrix2) <- c("All Autosomes", "all saDMPs", "Hypermethylated in females", "Hypermethylated in males")

results_matrix2[1,1] <- sum(backgroundFeature$Freq[1])
results_matrix2[1,2] <- sum(backgroundFeature$Freq[2] + sum(backgroundFeature$Freq[5]))
results_matrix2[1,3] <- sum(backgroundFeature$Freq[4])
results_matrix2[1,4] <- sum(backgroundFeature$Freq[3] + sum(backgroundFeature$Freq[6]))

results_matrix2[2,1] <- sum(saDMPFeature$Freq[1])
results_matrix2[2,2] <- sum(saDMPFeature$Freq[2] + sum(saDMPFeature$Freq[5]))
results_matrix2[2,3] <- sum(saDMPFeature$Freq[4])
results_matrix2[2,4] <- sum(saDMPFeature$Freq[3] + sum(saDMPFeature$Freq[6]))

results_matrix2[3,1] <- sum(femaleFeature$Freq[1])
results_matrix2[3,2] <- sum(femaleFeature$Freq[2] + sum(femaleFeature$Freq[5]))
results_matrix2[3,3] <- sum(femaleFeature$Freq[4])
results_matrix2[3,4] <- sum(femaleFeature$Freq[3] + sum(femaleFeature$Freq[6]))

results_matrix2[4,1] <- sum(maleFeature$Freq[1])
results_matrix2[4,2] <- sum(maleFeature$Freq[2] + sum(maleFeature$Freq[5]))
results_matrix2[4,3] <- sum(maleFeature$Freq[4])
results_matrix2[4,4] <- sum(maleFeature$Freq[3] + sum(maleFeature$Freq[6]))


hmap2 <- results_matrix2
results_matrix2 <- t(results_matrix2)
results_matrix_aves2 <- results_matrix2

results_matrix_aves2[,1] <- results_matrix2[,1]/sum(results_matrix2[,1])
results_matrix_aves2[,2] <- results_matrix2[,2]/sum(results_matrix2[,2])
results_matrix_aves2[,3] <- results_matrix2[,3]/sum(results_matrix2[,3])
results_matrix_aves2[,4] <- results_matrix2[,4]/sum(results_matrix2[,4])

rownames(results_matrix_aves2)<-c("CpG islands","Shelf","Open sea regions","Shore")
colnames(results_matrix_aves2)<-c("Autosomes","saDMPs","HM in females","HM in males")


pdf("/home/og16379/diff_cpg_fm/pdf/stackedBarplot2Rep.pdf")
par(mar=c(4,4,4,10.5))
par(xpd = NA)
barplot(results_matrix_aves2, col=cbPalette, yaxt = "none",cex.names=.4)
axis(2, seq(0, 1, 0.1), las=2,labels = paste0(seq(0, 100, 10), "%"))
legend(5,0.75,rownames(results_matrix_aves2),fill=cbPalette, bty = "n")
dev.off()


hmap_aves2 <- hmap2


hmap_aves2[1,] <- hmap2[1,]/sum(hmap2[1,])
hmap_aves2[2,] <- hmap2[2,]/sum(hmap2[2,])
hmap_aves2[3,] <- hmap2[3,]/sum(hmap2[3,])
hmap_aves2[4,] <- hmap2[4,]/sum(hmap2[4,])


hmap_comp2 <- hmap_aves2[2:4,]
hmap_comp2[1,] <- log2(hmap_comp2[1,]/hmap_aves2[1,])
hmap_comp2[2,] <- log2(hmap_comp2[2,]/hmap_aves2[1,])
hmap_comp2[3,] <- log2(hmap_comp2[3,]/hmap_aves2[1,])


hmap_comp[hmap_comp > 2] <- 2
hmap_comp[hmap_comp < -2] <- -2

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
custom_at <- seq(-2.5, 2.5, by = (5/30))
cols_contrast <- rev(c(colorRampPalette(c(cbbPalette[7], "white", cbbPalette[6]))(30)))

rownames(hmap_comp2)<-c("saDMPs","HM in females","HM in males")
colnames(hmap_comp2) <-c("CpG islands","Shelf","Open sea regions","Shore")
pdf("/home/og16379/diff_cpg_fm/pdf/heatmap2Rep.pdf",width=5)
  levelplot(t(hmap_comp2),
    at = custom_at,
    xlab = "log2 (obs/exp)",ylab="",
    scales=list(x=list(rot=90)),
    col.regions = cols_contrast)
dev.off()
