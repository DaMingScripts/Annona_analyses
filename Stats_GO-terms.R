setwd("~/Dropbox/HUB/4._Annona_Genome/RESULTS/OrthoMCL")
library(stats)
library(devtools)
library(ggord)


data<-read.delim(file="data_GO-term_occurences_for_stats.csv", header=TRUE, sep=",")
data$case = factor(data$case)
data2<-data.frame(data,row.names=data$case)
data2$case<-NULL


######## BP95
## Backgound vs Dicots_95BP
data_backgroundvsdicot95 <- data2[c(4, 7), ]
fisher.test(data_backgroundvsdicot95, workspace=2e+07,hybrid=TRUE, conf.int = TRUE, conf.level=0.95, simulate.p.value = TRUE)

## Backgound vs BasalAnnona_95BP
data_backgroundvsbasalAnnona95 <- data2[c(5, 7), ]
fisher.test(data_backgroundvsbasalAnnona95, workspace=2e+07,hybrid=TRUE, conf.int = TRUE, conf.level=0.95, simulate.p.value = TRUE)

## Backgound vs AnnMonocot_95BP
data_backgroundvsAnnMonocot95 <- data2[c(6, 7), ]
fisher.test(data_backgroundvsAnnMonocot95, workspace=2e+07,hybrid=TRUE, conf.int = TRUE, conf.level=0.95, simulate.p.value = TRUE)

## BasalAnnona_95BP vs Dicots_95BP
data_basalAnnona95vsdicot95 <- data2[c(4, 5), ]
fisher.test(data_basalAnnona95vsdicot95, workspace=2e+07,hybrid=TRUE, conf.int = TRUE, conf.level=0.95, simulate.p.value = TRUE)

## BasalAnnona_95BP vs AnnMonocot_BP95
data_AnnMonocot95vsbasalAnnona95 <- data2[c(5, 6), ]
fisher.test(data_AnnMonocot95vsbasalAnnona95, workspace=2e+07,hybrid=TRUE, conf.int = TRUE, conf.level=0.95, simulate.p.value = TRUE)



####### BP 70
## Backgound vs Dicots_70BP
data_backgroundvsdicot70 <- data2[c(3, 7), ]
fisher.test(data_backgroundvsdicot70, workspace=2e+07,hybrid=TRUE, conf.int = TRUE, conf.level=0.95, simulate.p.value = TRUE)

## Backgound vs BasalAnnona_70BP
data_backgroundvsbasalAnnona70 <- data2[c(2, 7), ]
fisher.test(data_backgroundvsbasalAnnona70, workspace=2e+07,hybrid=TRUE, conf.int = TRUE, conf.level=0.95, simulate.p.value = TRUE)

## Backgound vs AnnMonocot_70BP
data_backgroundvsAnnMonocot70 <- data2[c(1, 7), ]
fisher.test(data_backgroundvsAnnMonocot70, workspace=2e+07,hybrid=TRUE, conf.int = TRUE, conf.level=0.95, simulate.p.value = TRUE)

## BasalAnnona_70BP vs Dicots_70BP
data_basalAnnona70vsdicot70 <- data2[c(3, 5), ]
fisher.test(data_basalAnnona70vsdicot70, workspace=2e+07,hybrid=TRUE, conf.int = TRUE, conf.level=0.95, simulate.p.value = TRUE)

## BasalAnnona_70BP vs AnnMonocot_BP70
data_AnnMonocot70vsbasalAnnona70 <- data2[c(2, 1), ]
fisher.test(data_AnnMonocot70vsbasalAnnona70, workspace=2e+07,hybrid=TRUE, conf.int = TRUE, conf.level=0.95, simulate.p.value = TRUE)




PCA_data2<-prcomp(data2, center = TRUE, scale. = TRUE) 
summary(PCA_data2)
plot(PCA_data2, type="l")
data_names<-c("AnnMonocotBP70","basalAnnona_BP70","DicotsBP70","DicotsBP95","basalAnnona_BP95","AnnMonocotBP95","background")
data_names2<-c(1,2,3,4,5,6,7)
load_names<-c("antioxidant_activity","catalytic_activity")
pts_PCA<-PCA_data2$x
pts_PCA <- subset(pts_PCA, select = PC1:PC2)
pts_PCA <- as.data.frame(pts_PCA)
ggord(PCA_data2, obslab = TRUE, size = 5, facet=TRUE, nfac=2, addpts = pts_PCA, repel = TRUE)

###########################






vsd <- varianceStabilizingTransformation(dds)
data <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
plotPCA <- ggplot(data, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=names),hjust=0.25, vjust=-0.5, show_guide = F)
ggsave("PCA.pdf", plot = plotPCA)