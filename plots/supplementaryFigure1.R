################################################################
# Characterising Autosomal sex differences between males and females
# Olivia Grant
###############################################################

# below is all the code to produce the figures shown in Supplementary figure 1A-1D
# previous code must have been run in order to produce these figures

###############################################################
# first work with discovery data set
###############################################################

pvals<- dmp$P.Value

 dmpQQ <- function(pvals){

            fold.range <- c(0,6)
             prob.range <- c(0,300)
             plot(
                     x = 0, y = 0, xlim = fold.range, ylim = prob.range,
                     xlab = "", ylab = "",
                     axes = TRUE, type = "n", cex.axis = 1.0, cex.lab = 1.0)

             coords <- par("usr")
             gx <- grconvertX(coords[1:2], "user", "in")
             gy <- grconvertY(coords[3:4], "user", "in")
             width <- max(gx) - min(gx)
             height <- max(gy) - min(gy)

             tmp <- tempfile()
             png(tmp, width = width, height = height, units = "in", res = 500, bg = "transparent")

             plot.new()
             par(mar = c(0.0, 0.0, 0.0, 0.0))
             qq(pvals)
             dev.off()
             panel <- readPNG(tmp)
             rasterImage(panel, coords[1], coords[3], coords[2], coords[4])
             mtext(text = "Expected -log10 (p)", side = 1, line = 4, cex = 1)
             mtext(text = "Observed-log10(p)", side = 2, line = 2.5, cex = 1)

    }

   pdf("/home/og16379/diff_cpg_fm/pdf/rasterisedQQ.pdf",width=8)
     dmpQQ(pvals)
 dev.off()

 P_lambda(pvals)


 ################################################################################
 # Generate cell type and perform wilcoxon test
 ################################################################################

 #first estimate cell type propotions
 ecc<-estimateCellCounts.gds(gfile)
 #remove outliers
 eccDf<-as.data.frame(ecc)
 eccDf<-eccDf[info$barcode,]
 eccDf$Sex<-info$Sex

  x<- which(eccDf$Sex=="")
  eccDf <- eccDf[-x]

 eccDf$Sex<- gsub('F', 'Females', eccDf$Sex)
 eccDf$Sex<- gsub('M', 'Males', eccDf$Sex)

 eccMelt<-melt(eccDf)

 #plot
 ggplot(eccMelt,aes(x=variable,y=value,fill=Sex))+
 geom_boxplot()+
 labs(x="Blood cell types",y="Proportion")+theme(legend.title=element_blank())+
   theme(axis.line = element_line(colour = "black"),
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     panel.background=element_rect(fill="white",colour="white"),
     panel.border = element_rect(fill=NA,color="black"))+
     scale_fill_manual(values=c(femaleColour,maleColour))

 #wilcoxon p values
 f<-split(eccMelt,eccMelt$Sex)
 femaleData<-f$Females
 maleData<-f$Males
 maleData<-split(maleData,maleData$variable)
 femaleData<-split(femaleData,femaleData$variable)

 maleMatrix<-matrix(0,1132,6)
 colnames(maleMatrix)<-c("CD8T","CD4T","NK","Bcell","Mono","Gran")
 maleMatrix[,1]<-maleData$CD8T$value
 maleMatrix[,2]<-maleData$CD4T$value
 maleMatrix[,3]<-maleData$NK$value
 maleMatrix[,4]<-maleData$Bcell$value
 maleMatrix[,5]<-maleData$Mono$value
 maleMatrix[,6]<-maleData$Gran$value

 femaleMatrix<-matrix(0,1339,6)
 colnames(femaleMatrix)<-c("CD8T","CD4T","NK","Bcell","Mono","Gran")
 femaleMatrix[,1]<-femaleData$CD8T$value
 femaleMatrix[,2]<-femaleData$CD4T$value
 femaleMatrix[,3]<-femaleData$NK$value
 femaleMatrix[,4]<-femaleData$Bcell$value
 femaleMatrix[,5]<-femaleData$Mono$value
 femaleMatrix[,6]<-femaleData$Gran$value

 #repeat in loop
 wilcox.test(femaleMatrix[,1],maleMatrix[,1])



 
 ###############################################################
 # first work with validation data set
 ###############################################################

 pvals2<- dmp2$P.Value

  dmpQQ2 <- function(pvals2){

             fold.range <- c(0,6)
              prob.range <- c(0,300)
              plot(
                      x = 0, y = 0, xlim = fold.range, ylim = prob.range,
                      xlab = "", ylab = "",
                      axes = TRUE, type = "n", cex.axis = 1.0, cex.lab = 1.0)

              coords <- par("usr")
              gx <- grconvertX(coords[1:2], "user", "in")
              gy <- grconvertY(coords[3:4], "user", "in")
              width <- max(gx) - min(gx)
              height <- max(gy) - min(gy)

              tmp <- tempfile()
              png(tmp, width = width, height = height, units = "in", res = 500, bg = "transparent")

              plot.new()
              par(mar = c(0.0, 0.0, 0.0, 0.0))
              qq(pvals2)
              dev.off()
              panel <- readPNG(tmp)
              rasterImage(panel, coords[1], coords[3], coords[2], coords[4])
              mtext(text = "Expected -log10 (p)", side = 1, line = 4, cex = 1)
              mtext(text = "Observed-log10(p)", side = 2, line = 2.5, cex = 1)

     }

    pdf("/home/og16379/diff_cpg_fm/pdf/rasterisedQQ2.pdf",width=8)
      dmpQQ2(pvals2)
  dev.off()

  P_lambda(pvals2)


  ################################################################################
  # Generate cell type and perform wilcoxon test
  ################################################################################

  #first estimate cell type propotions
  ecc2<-estimateCellCounts.gds(gfile_2)
  #remove outliers
  eccDf2<-as.data.frame(ecc2)
  eccDf2<-eccDf2[info2$barcode,]
  eccDf2$Sex<-info2$Sex

   x<- which(eccDf2$Sex=="")
   eccDf2 <- eccDf2[-x]

  eccDf2$Sex<- gsub('F', 'Females', eccDf2$Sex)
  eccDf2$Sex<- gsub('M', 'Males', eccDf2$Sex)

  eccMelt2<-melt(eccDf2)

  #plot
  ggplot(eccMelt2,aes(x=variable,y=value,fill=Sex))+
  geom_boxplot()+
  labs(x="Blood cell types",y="Proportion")+theme(legend.title=element_blank())+
    theme(axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background=element_rect(fill="white",colour="white"),
      panel.border = element_rect(fill=NA,color="black"))+
      scale_fill_manual(values=c(femaleColour,maleColour))

  #wilcoxon p values
  f2<-split(eccMelt2eccMelt2$Sex)
  femaleData2<-f2$Females
  maleData2<-f2$Males
  maleData2<-split(maleData2,maleData2$variable)
  femaleData2<-split(femaleData2,femaleData2$variable)

  maleMatrix2<-matrix(0,1132,6)
  colnames(maleMatrix2)<-c("CD8T","CD4T","NK","Bcell","Mono","Gran")
  maleMatrix2[,1]<-maleData2$CD8T$value
  maleMatrix2[,2]<-maleData2$CD4T$value
  maleMatrix2[,3]<-maleData2$NK$value
  maleMatrix2[,4]<-maleData2$Bcell$value
  maleMatrix2[,5]<-maleData2$Mono$value
  maleMatrix2[,6]<-maleData2$Gran$value

  femaleMatrix2<-matrix(0,1339,6)
  colnames(femaleMatrix2)<-c("CD8T","CD4T","NK","Bcell","Mono","Gran")
  femaleMatrix2[,1]<-femaleData2$CD8T$value
  femaleMatrix2[,2]<-femaleData2$CD4T$value
  femaleMatrix2[,3]<-femaleData2$NK$value
  femaleMatrix2[,4]<-femaleData2$Bcell$value
  femaleMatrix2[,5]<-femaleData2$Mono$value
  femaleMatrix2[,6]<-femaleData2$Gran$value

  #repeat in loop
  wilcox.test(femaleMatrix2[,1],maleMatrix2[,1])

  ################################################################################
  # Generate venn diagram 
  ################################################################################

  
  if (!require(devtools)) install.packages("devtools")
  devtools::install_github("yanlinlin82/ggvenn")

  library(ggvenn)

  ################################################################################
  # Generate Figure 2A (venn)
  ################################################################################

  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  maleColour<-cbbPalette[7]
  femaleColour<-cbbPalette[6]

  x <- list(
    "Discovery saDMPs" = dasenAutosomeCleaned$Row.names,
    "Validation saDMPs" = dasenAutosomeCleaned2$Row.names,
    "Significant Discovery saDMPs" = allSigProbes$Row.names,
    "Significant Validation saDMPs" = allSigProbes2$Row.names
    )

  pdf("/home/og16379/diff_cpg_fm/pdf/venn.pdf")
  ggvenn(
    x,
    fill_color = cbPalette,
  show_percentage=FALSE,
  stroke_color="white",
    stroke_size = 0.5, set_name_size = 1.8
    )
  dev.off()
