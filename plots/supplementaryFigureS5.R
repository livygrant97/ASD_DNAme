################################################################
# Characterising Autosomal sex differences between males and females
# Olivia Grant
###############################################################

# below is all the code to produce the figures shown in Supplementary figure 5A-5C
# previous code must have been run in order to produce these figures

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biomaRt")

##################################################
# 3 volcano plots #
# one for sex chr, one for autosomes, one for sadmps #
##################################################
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
IDs <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
                 filters = "ensembl_gene_id",values=rownames(res),
                 mart = mart)


saDMPIds <- getBM(attributes = c("external_gene_name","ensembl_gene_id"),
                 filters = "external_gene_name",values=replicatedResults$gene,
                 mart = mart)


notMapped<-rownames(res)
notMapped<-notMapped[!notMapped %in% IDs$ensembl_gene_id]
notMapped<-data.frame("ensembl_gene_id"=notMapped,"external_gene_name"="")

geneInfo<-rbind(IDs,notMapped)
geneInfo1<-geneInfo[match(rownames(res), geneInfo$ensembl_gene_id), ]
res$genename <- geneInfo1$external_gene_name

resDf <- as.data.frame(res)
sigRes <- subset(resDf,padj <= 0.05)

#first, subset count matrix according to groups above
countMatrixSexChr<-subset(countMatrix, Chromosome.Name == "chrX" | Chromosome.Name =="chrY")
countMatrixAutosomes<- countMatrix %>% dplyr::filter(Chromosome.Name != "chrX" & Chromosome.Name !="chrY")
countMatrixDMPs<-subset(countMatrix, Geneid %in% replicatedResults$gene)

#subset results according to above
resultsSexChr<-subset(res,rownames(res) %in% countMatrixSexChr$Ensembl.Gene.ID)
resultsAutosomes<-subset(res,rownames(res) %in% countMatrixAutosomes$Ensembl.Gene.ID)
resultsDMPs<-subset(res,rownames(res) %in% countMatrixDMPs$Ensembl.Gene.ID)



#make volcano plot for each
pdf("/home/og16379/diff_cpg_fm/pdf/VolcanoSexChr.pdf",width=10)
v1<-EnhancedVolcano(resultsSexChr,
  lab = resultsSexChr$genename,
  selectLab = resultsSexChr$genename[1:30],
  x = 'log2FoldChange',
  y = 'padj',
  xlab = bquote(~Log[2]~ 'fold change'),
  FCcutoff = 1.0,
  colAlpha = 1,
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 4.0,
  widthConnectors = 0.75)
  v1
dev.off()

pdf("/home/og16379/diff_cpg_fm/pdf/VolcanoAutosomes.pdf",width=10)
v2<-EnhancedVolcano(resultsAutosomes,
  lab=resultsAutosomes$genename,
  selectLab = resultsAutosomes$genename[1:30],
  x = 'log2FoldChange',
  y = 'padj',
  xlab = bquote(~Log[2]~ 'fold change'),
  FCcutoff = 1.0,
  colAlpha = 1,
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 4.0,
  widthConnectors = 0.2,)
v2
dev.off()

pdf("/home/og16379/diff_cpg_fm/pdf/VolcanoSaDMPs.pdf",width=10)
v3<-EnhancedVolcano(resultsDMPs,
  lab=resultsDMPs$genename,
  xlim=c(-1,1.7),
  selectLab = resultsDMPs$genename[1:50],
  x = 'log2FoldChange',
  y = 'padj',
  xlab = bquote(~Log[2]~ 'fold change'),
  colAlpha = 1,
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 4.0,
  widthConnectors = 0.2,)
v3
dev.off()
