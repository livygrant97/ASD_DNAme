################################################################
# Characterising Autosomal sex differences between males and females
# Olivia Grant
###############################################################

#below is the code to produce all figures in figure 2

# load packages

if(!require("VennDiagram", character.only = TRUE)){
  install.packages("VennDiagram")
}
library(VennDiagram)


if(!require("seqLogo", character.only = TRUE)){
  install.packages("seqLogo")
}
library(seqLogo)

if(!require("clusterProfiler", character.only = TRUE)){
  install.packages("clusterProfiler")
}
library(clusterProfiler)


#################################################################################
#### Analysing the data ####
#################################################################################
#colours we use for plotting throughout the manuscript
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
maleColour<-cbbPalette[7]
femaleColour<-cbbPalette[6]


################################################################################
# Generate Figure 2A (venn)
################################################################################

pdf("/home/og16379/diff_cpg_fm/pdf/overlappingTFRep.pdf")
grid.newpage()
draw.pairwise.venn(area1=length(enrichedInDMPsHMInFemales),area2=length(enrichedInDMPsHMInMales),cross.area=length(overlap),fill=c(femaleColour,maleColour),lty="blank",category=c("Female","Male"),cex=1)
dev.off()

################################################################################
# Generate Figure 2B (motifs)
################################################################################

#query databases
SOX9MDB<-query(query(query(MotifDb,"hsapiens"), "SOX9"), "HOCOMOCOv10")
SRYMDB<-query(query(query(MotifDb,"hsapiens"), "SRY"), "HOCOMOCOv10")
ESR1MDB<-query(query(query(MotifDb,"hsapiens"), "ESR1"), "HOCOMOCOv10")
#get motifs
SOX9MDB<-SOX9MDB[[1]]
SRY9MDB<-SRYMDB[[1]]
ESR1MDB<-ESR1MDB[[1]]

#SRY
pdf("/home/og16379/diff_cpg_fm/pdf/SRYLogo.pdf", height=5)
seqLogo(SRY9MDB)
dev.off()

#SOX9
pdf("/home/og16379/diff_cpg_fm/pdf/SOX9Logo.pdf", height=5)
seqLogo(SOX9MDB)
dev.off()

#ESR1
pdf("/home/og16379/diff_cpg_fm/pdf/ESR1Logo.pdf", height=5)
seqLogo(ESR1MDB)
dev.off()

################################################################################
# next, perform GO and KEGG enrichment on enriched motifs
################################################################################

#convert to entrez ids
male.df <- bitr(maleOnly, fromType = "SYMBOL",
        toType = c("ENTREZID", "SYMBOL","ENSEMBL"),
        OrgDb = org.Hs.eg.db)
head(male.df)

#convert to entrez ids
female.df <- bitr(femaleOnly, fromType = "SYMBOL",
        toType = c("ENTREZID", "SYMBOL","ENSEMBL"),
        OrgDb = org.Hs.eg.db)
head(female.df)

#perform kegg enrichment
maleKegg<-enrichKEGG(male.df$ENTREZID)
femaleKegg<-enrichKEGG(female.df$ENTREZID)

#perform go enrichment
maleGO<-enrichGO(male.df$ENTREZID,org.Hs.eg.db)
femaleGO<-enrichGO(female.df$ENTREZID,org.Hs.eg.db)

################################################################################
# Generater Figure 2C
################################################################################

#the only significant KEGG output was for females
pdf("/home/og16379/diff_cpg_fm/pdf/femaleKeggEnrichmentRep.pdf", height=5)
dotplot(femaleKegg,showCategory=25)
dev.off()

################################################################################
# Generate Supplementary Figure S3
################################################################################

#S2A

pdf("/home/og16379/diff_cpg_fm/pdf/femaleGOenrichmentRep.pdf",width=20)
dotplot(femaleGO,showCategory=27)
dev.off()

#S2B
pdf("/home/og16379/diff_cpg_fm/pdf/maleGOenrichmentRep.pdf",width=20)
dotplot(maleGO,showCategory=16)
dev.off()


################################################################################
# figure 2D and 2E were produced in cytoscape
# https://cytoscape.org/
################################################################################
