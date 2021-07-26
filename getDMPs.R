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
#get methylated and unmethylated
mns <- methylated(gfile)[,]
uns <- unmethylated(gfile)[,]
#remove sex chromosome probes
autosomes <- !(fData(gfile)$CHR %in% c('X', 'Y'))
print("Dasen normalization on CpGs from autosomes...")
dasen_autosome <- dasen(mns[autosomes, ], uns[autosomes, ], fData(gfile)$Infinium_Design_Type[autosomes])
#write to file
fwrite(data.frame(ID_REF=rownames(dasen_autosome), dasen_autosome, check.names=F), file=save_file, sep='\t')suppressMessages(require(bigmelon))


# next, we can calculate dmps using champ dmp & t test method
# load normalised betas
suppressMessages(require(ChAMP))   ## champ.DMP
suppressMessages(require(data.table))
#arg <- commandArgs(T)
autosome <- TRUE
#autosome <- arg[1]
gds_file <- '/storage/st05d/Exeter/UnderSocMeth/USM1/USM_WF.gds'
save_file <- 'diff_cpg_between_F_M.xls'
gfile <- openfn.gds(gds_file)
if(autosome){
    print("Find differential methylated CpG sites on autosomes")
    ## load normalized autosomal beta values
    normed_beta <- '/path/to/file/autosome_norm_beta.xls'
    # betas <- index.gdsn(gfile, 'autobetas')[1:847214,]
    betas <- fread(normed_beta, data.table=FALSE)
    rownames(betas) <- betas$ID_REF
    betas <- betas[, -1]
} else {
    betas <- betas(gfile)[,]
}

## get sex annotation
info <- pData(gfile)[, c('barcode', 'nsex')]
info$nsex <- gsub('2', 'F', info$nsex)
info$nsex <- gsub('1', 'M', info$nsex)

# remove samples identified as outliers using the DNAme based sex classifer (Wang et al 2020)
# https://www.biorxiv.org/content/10.1101/2020.10.19.345090v1

outlier <- c('200611820013_R08C01', '200603220075_R08C01', '200864580018_R02C01', '200611820020_R07C01')
info <- info[!(info$barcode %in% outlier), ]
betas <- betas[, info$barcode]
pheno_group <- info$nsex


## Find Differential Methylation Positions  by F-test from limma
## add adj.Pval argument and set to 1 if you want all stats for all probes to make the volcano plot in figure 1
f_DMP <- champ.DMP(beta=betas, pheno=pheno_group, arraytype='EPIC')[[1]]

## Find Differential Methylation Positions  by T-test
print("Performing T test: ")
t_Pvals <- apply(betas, 1, function(x){t.test(as.numeric(x)~pheno_group)$p.value})  ## extremly time consuming
t_Padj <- p.adjust(t_Pvals, method = 'BH')
t_DMP <- data.frame(t_Pvals, t_Padj)
t_DMP <- t_DMP[t_DMP$t_Padj < 0.05, ]

#merge the two results together
dmp <- merge(f_DMP, t_DMP, by=0, sort=FALSE)
fwrite(dmp, file=save_file, sep='\t')



################################################################
# NOW, working with the validation data set
###############################################################

#first, we must normalise our data with the autosomes only
#load packages
suppressMessages(require(bigmelon))
suppressMessages(require(data.table))

#file name
save_file <- 'autosome_norm_beta2.xls'
#path to gds file or EPIC file
gds_file2 <- '/storage/st05d/Exeter/UnderSocMeth/ISER/gds_files_corrected/UKHLS_2_NormalisedBetas.gds'
#open gds file
gfile_2 <- openfn.gds(gds_file2)
#get methylated and unmethylated
#mns <- methylated(gfile_2)[,]
#uns <- unmethylated(gfile_2)[,]
betas <- betas(gfile_2)
#remove sex chromosome probes
fdata<-as.data.frame(read.gdsn(index.gdsn(gfile_2,"fData")))
autosomes <- !(fdata$CHR %in% c('chrX', 'chrY'))
print("Dasen normalization on CpGs from autosomes...")
dasen_autosome_2 <- dasen(betas[autosomes], fdata$Infinium_Design_Type[autosomes])
save(dasen_autosome2,file="/home/og16379/diff_cpg_fm/data/dasen_autosome2.Rdata")
#write to file
fwrite(data.frame(ID_REF=rownames(dasen_autosome2), dasen_autosome2, check.names=F), file=save_file, sep='\t')


# next, we can calculate dmps using champ dmp & t test method
# load normalised betas
suppressMessages(require(ChAMP))   ## champ.DMP
suppressMessages(require(data.table))
suppressMessages(require(bigmelon))
#arg <- commandArgs(T)
autosome <- TRUE
#autosome <- arg[1]

## get sex annotation
info <- pData(gfile_2)[, c('Sample','Rack.Barcode', 'Sex')]
info$Sex <- gsub('2', 'F', info$Sex)
info$Sex <- gsub('1', 'M', info$Sex)

# remove samples identified as outliers using the DNAme based sex classifer (Wang et al 2020)
# https://www.biorxiv.org/content/10.1101/2020.10.19.345090v1

outlier <- c('203991410050_R06C01', '203991410050_R07C01', '204026590012_R07C01', '203998240051_R07C01',
    '204026590078_R04C01', '203994670097_R08C01', '203960330159_R08C01', '204022160070_R03C01',
    '203994670097_R02C01', '203991470090_R05C01', '204026590012_R06C01')
    
info <- info[!(info$barcode %in% outlier), ]
dasen_autosome2 <- dasen_autosome2[, info$barcode]
pheno_group <- info$Sex


## Find Differential Methylation Positions  by F-test from limma
## add adj.Pval argument and set to 1 if you want all stats for all probes to make the volcano plot in figure 1
f_DMP <- champ.DMP(beta=dasen_autosome2, pheno=pheno_group, adjPVal=1,arraytype='EPIC')[[1]]

## Find Differential Methylation Positions  by T-test
print("Performing T test: ")
t_Pvals <- apply(dasen_autosome2, 1, function(x){t.test(as.numeric(x)~pheno_group)$p.value})  ## extremly time consuming
t_Padj <- p.adjust(t_Pvals, method = 'BH')
t_DMP <- data.frame(t_Pvals, t_Padj)
t_DMP <- t_DMP[t_DMP$t_Padj < 0.05, ]

#merge the two results together
dmp2 <- merge(f_DMP, t_DMP, by=0, sort=FALSE)
fwrite(dmp2, file=save_file, sep='\t')

