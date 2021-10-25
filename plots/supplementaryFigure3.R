################################################################################
# Generate Supplementary Figure S3
################################################################################

#S3A

pdf("/home/og16379/diff_cpg_fm/pdf/femaleGOenrichment.pdf")
dotplot(femaleGO,showCategory=50)
dev.off()

#S3B
pdf("/home/og16379/diff_cpg_fm/pdf/maleGOenrichment.pdf")
dotplot(maleGO,showCategory=50)
dev.off()
