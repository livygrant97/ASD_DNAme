################################################################
# Characterising Autosomal sex differences between males and females
# Olivia Grant
###############################################################

################################################################
# FIRST, working with the discovery data set
###############################################################

#first, we must normalise our data with the autosomes only
#load packages
suppressMessages(require(bigmelon))
suppressMessages(require(data.table))

#file name
save_file <- 'autosome_norm_beta.xls'
#path to gds file or EPIC file
gds_file <- '/storage/projects/Exeter/UnderSocMeth/USM1/USM_WF.gds'
#open gds file
gfile <- openfn.gds(gds_file)
#get methylated and unmethylated and betas
mns <- methylated(gfile)[,]
uns <- unmethylated(gfile)[,]
betas <- betas(gfile)
#remove sex chromosome probes
autosomes <- !(fData(gfile)$CHR %in% c('X', 'Y'))
# perform dasen normalisation on autosomes only to avoid
# introducing technical bias from sex chromosomes
print("Dasen normalization on CpGs from autosomes...")
dasen_autosome <- dasen(mns[autosomes, ], uns[autosomes, ], fData(gfile)$Infinium_Design_Type[autosomes])
# save normalised matrix
save(dasen_autosome,file="/home/og16379/diff_cpg_fm/data/firstDatasetNormalised.RData")

# next, we can calculate dmps using champ dmp & t test method
# load packages required
suppressMessages(require(ChAMP))   ## champ.DMP
suppressMessages(require(data.table))

# get sex annotation
info <- pData(gfile)[, c('barcode', 'nsex')]
# change class
info$nsex <- gsub('2', 'F', info$nsex)
info$nsex <- gsub('1', 'M', info$nsex)

# remove samples identified as outliers using the DNAme based sex classifer (Wang et al 2020)
# https://www.biorxiv.org/content/10.1101/2020.10.19.345090v1
# assign vector of outlier barcodes
outlier <- c('200611820013_R08C01', '200603220075_R08C01', '200864580018_R02C01', '200611820020_R07C01')
# remove misclassified barcodes from data
info <- info[!(info$barcode %in% outlier), ]
dasen_autosome <- dasen_autosome[, info$barcode]
pheno_group <- info$nsex

## Find Differential Methylation Positions  by F-test from limma
## add adj.Pval argument and set to 1 if you want all stats for all probes to make the volcano plot in figure 1
f_DMP <- champ.DMP(beta=dasen_autosome, pheno=pheno_group, arraytype='EPIC')[[1]]
## Find Differential Methylation Positions  by T-test
print("Performing T test: ")
t_Pvals <- apply(dasen_autosome, 1, function(x){t.test(as.numeric(x)~pheno_group)$p.value})  ## extremly time consuming
t_Padj <- p.adjust(t_Pvals, method = 'BH')
t_DMP <- data.frame(t_Pvals, t_Padj)
t_DMP <- t_DMP[t_DMP$t_Padj < 0.05, ]

#merge the two results together
dmp <- merge(f_DMP, t_DMP, by=0, sort=FALSE)
#save the results
save(dmp,file="/home/og16379/diff_cpg_fm/data/dmpFirstRound.RData")


################################################################
# NOW, working with the validation data set
###############################################################
#first, we must normalise our data with the autosomes only again
#load packages
suppressMessages(require(bigmelon))
suppressMessages(require(data.table))

#path to gds file or EPIC file
gds_file2 <- '/storage/st05d/Exeter/UnderSocMeth/ISER/USM2_WF.gds'
#open gds file
gfile_2 <- openfn.gds(gds_file2)
#get methylated and unmethylated
mns2 <- methylated(gfile_2)[,]
uns2 <- unmethylated(gfile_2)[,]
# get f data
fdata<-as.data.frame(read.gdsn(index.gdsn(gfile_2,"fData")))
autosomes2 <- !(fdata$chr %in% c('chrX', 'chrY'))
print("Dasen normalization on CpGs from autosomes...")
# perform dasen normalisation on autosomes only to avoid
# introducing technical bias to autosomal cpgs
dasen_autosome2 <- dasen(mns2[autosomes2, ], uns2[autosomes2, ], fdata$Type[autosomes2])
# save to rdata
save(dasen_autosome2,file="/home/og16379/diff_cpg_fm/data/dasen_autosome2round.Rdata")

# next, we can calculate dmps using champ dmp & t test method
# load packages
suppressMessages(require(ChAMP))   ## champ.DMP
suppressMessages(require(data.table))
suppressMessages(require(bigmelon))

## get sex annotation
info2 <- pData(gfile_2)[, c('barcode','Sex')]
info2$Sex <- gsub('2', 'F', info2$Sex)
info2$Sex <- gsub('1', 'M', info2$Sex)

# remove samples identified as outliers using the DNAme based sex classifer (Wang et al 2020)
# https://www.biorxiv.org/content/10.1101/2020.10.19.345090v1

outlier2 <- c('203991410050_R06C01', '203991410050_R07C01', '204026590012_R07C01', '203998240051_R07C01',
    '204026590078_R04C01', '203994670097_R08C01', '203960330159_R08C01', '204022160070_R03C01',
    '203994670097_R02C01', '203991470090_R05C01', '204026590012_R06C01')

info2 <- info2[!(info2$barcode %in% outlier2), ]
dasen_autosome2 <- dasen_autosome2[, info2$barcode]
pheno_group2 <- info2$Sex


## Find Differential Methylation Positions  by F-test from limma
## add adj.Pval argument and set to 1 if you want all stats for all probes to make the volcano plot in figure 1
f_DMP2 <- champ.DMP(beta=dasen_autosome2, pheno=pheno_group2,arraytype='EPIC')[[1]]

## Find Differential Methylation Positions  by T-test
print("Performing T test: ")
t_Pvals2 <- apply(dasen_autosome2, 1, function(x){t.test(as.numeric(x)~pheno_group2)$p.value})  ## extremly time consuming
t_Padj2 <- p.adjust(t_Pvals2, method = 'BH')
t_DMP2 <- data.frame(t_Pvals2, t_Padj2)
t_DMP2 <- t_DMP2[t_DMP2$t_Padj2 < 0.05, ]

#merge the two results together
dmp2 <- merge(f_DMP2, t_DMP2, by=0, sort=FALSE)
save(dmp2,file="/home/og16379/diff_cpg_fm/data/dmpSecondRound.RData")
fwrite(dmp2, file=save_file, sep='\t')
