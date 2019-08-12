setwd("~/Dropbox/HUB/4._Annona_Genome/RESULTS/OrthoMCL")
library(stats)
library(gplots)
library(ggord)
library(ggplot2)
library (gridExtra)
library(grid)
require(graphics); require(grDevices)

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}




data<-read.delim(file="genes_count_per_scaffold_all_topology.csv", header=TRUE, sep=",")
row.names(data)<-data$scaffold_ID

data<-as.matrix(data)
data_heatmap_nbgene<-subset(data, select= nb_gene_background:nb_gene_Dicots)
data_heatmap_percnbgene<-subset(data, select= scaffold_length:perc_nb_genes_Dicots)
data_heatmap_densitygene<-subset(data, select= density_gene_background:density_gene_Dicots)
data_heatmap_percgene<-subset(data, select= perc_gene_background:perc_gene_Dicots)


class(data_heatmap_nbgene)<- "numeric"
class(data_heatmap_densitygene)<- "numeric"
class(data_heatmap_percgene)<- "numeric"

data_heatmap_diffgene<-subset(data, select= diff_gene_background:diff_gene_Dicots)
class(data_heatmap_diffgene)<- "numeric"

pdf(file="heatmaps_genetic_structure_topology_hypotheses_scafnames.pdf", width=12, height=8)

# data_heatmap[1:] <-as.numeric(data_heatmap[1:])

heatmap.2(data_heatmap_nbgene, scale="row", na.rm = TRUE, main = "Nb of genes per scaffold", xlab = "Hypotheses", ylab = "scaffolds", srtCol = 0, adjCol=c(0.7,0.5), cexCol = 1, labCol = c("Background", "AnnMonocots", "BasalAnnona", "Dicots"), labRow = "")
heatmap.2(data_heatmap_densitygene, scale="row", na.rm = TRUE, main = "density of genes per scaffold", xlab = "Hypotheses", ylab = "scaffolds", srtCol = 0, adjCol=c(0.7,0.5), cexCol = 1, labCol = c("Background", "AnnMonocots", "BasalAnnona", "Dicots"), labRow = "")
heatmap.2(data_heatmap_percgene, scale="row", na.rm = TRUE, main = "percentage of genes per scaffold", xlab = "Hypotheses", ylab = "scaffolds", srtCol = 0, adjCol=c(0.7,0.5), cexCol = 1, labCol = c("Background", "AnnMonocots", "BasalAnnona", "Dicots"), labRow = "")
heatmap.2(data_heatmap_diffgene, scale="column", na.rm = TRUE, main = "difference of genes density with the background - per scaffold", xlab = "Hypotheses", ylab = "scaffolds", srtCol = 0, adjCol=c(0.7,0.5), cexCol = 1, labCol = c("Background", "AnnMonocots", "BasalAnnona", "Dicots"), labRow = "")

dev.off()

hm_nbgenes <- heatmap.2(data_heatmap_nbgene, scale="row", na.rm = TRUE, main = "Nb of genes per scaffold", xlab = "Hypotheses", ylab = "scaffolds", srtCol = 0, adjCol=c(0.7,0.5), cexCol = 1, labCol = c("Background", "AnnMonocots", "BasalAnnona", "Dicots"), labRow = "")
hm_densitygenes <- heatmap.2(data_heatmap_densitygene, scale="row", na.rm = TRUE, main = "density of genes per scaffold", xlab = "Hypotheses", ylab = "scaffolds", srtCol = 0, adjCol=c(0.7,0.5), cexCol = 1, labCol = c("Background", "AnnMonocots", "BasalAnnona", "Dicots"), labRow = "")
hm_percgene <- heatmap.2(data_heatmap_percgene, scale="row", na.rm = TRUE, main = "percentage of genes per scaffold", xlab = "Hypotheses", ylab = "scaffolds", srtCol = 0, adjCol=c(0.7,0.5), cexCol = 1, labCol = c("Background", "AnnMonocots", "BasalAnnona", "Dicots"), labRow = "")
hm_diffgene <- heatmap.2(data_heatmap_diffgene, scale="column", na.rm = TRUE, main = "difference of genes density with the background - per scaffold", xlab = "Hypotheses", ylab = "scaffolds", srtCol = 0, adjCol=c(0.7,0.5), cexCol = 1, labCol = c("Background", "AnnMonocots", "BasalAnnona", "Dicots"), labRow = "")

nbgenes_value <- hm_nbgenes$carpet
densitygenes_value <- hm_densitygenes$carpet
percgene_value <- hm_percgene$carpet
diffgene_value <- hm_diffgene$carpet

write.table(nbgenes_value, file="heatmap_nbgenes.txt", row.names=TRUE, col.names=TRUE)
write.table(densitygenes_value, file="heatmap_densitygenes.txt", row.names=TRUE, col.names=TRUE)
write.table(percgene_value, file="heatmap_percgene.txt", row.names=TRUE, col.names=TRUE)
write.table(diffgene_value, file="heatmap_diffgene.txt", row.names=TRUE, col.names=TRUE)

### Scatter-plot of data


## Density graph
data_ggplot<-read.delim(file="genes_count_per_scaffold_all_topology_ggplot.csv", header=TRUE, sep=",")

# inset
a <- ggplot(data_ggplot, aes(x=density_gene_background, y=y_density, colour=data_ggplot[,6])) + geom_point() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(),plot.background = element_rect(fill = "white"), panel.background = element_rect(fill = "lightgrey", colour = "white", size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.05, linetype = 'solid', colour = "black"),legend.position = "none")
# large graph
b <- ggplot(data_ggplot, aes(x=density_gene_background, y=y_density, colour=data_ggplot[,6])) + geom_point() + xlab("density in background") + ylab("density of genes supporting a given topology") + geom_rug(position="jitter", size=.2) + xlim(0,2E-6) + ylim(0,2E-6) + theme(legend.position = c(0,1.03), legend.title=element_blank(), legend.justification = c(-0.1,1)) + theme(legend.key=element_blank())

grid.newpage()
vpb_ <- viewport(width = 1, height = 1, x = 0.5, y = 0.5)  # the larger map
vpa_ <- viewport(width = 0.4, height = 0.3, x = 0.8, y = 0.847)  # the inset in upper right
print(b, vp = vpb_)
print(a, vp = vpa_)








## Percentage of genes supporting a topology 
data_heatmap_percnbgene<-read.delim(file="perc_nb_genes.csv", header=TRUE, sep=",")

data_heatmap_percnbgene$scaffold_length <- data_heatmap_percnbgene$scaffold_length / 1000


levels(data_heatmap_percnbgene$hypothesis)[levels(data_heatmap_percnbgene$hypothesis)=="perc_nb_genes_AnnMonocots"] <- "AnnMonocots" 
levels(data_heatmap_percnbgene$hypothesis)[levels(data_heatmap_percnbgene$hypothesis)=="perc_nb_genes_basalAnnona"] <- "basalAnn" 
levels(data_heatmap_percnbgene$hypothesis)[levels(data_heatmap_percnbgene$hypothesis)=="perc_nb_genes_Dicots"] <- "Dicots" 


test_AnnMonocots<-subset(data_heatmap_percnbgene[ which(data_heatmap_percnbgene$hypothesis=='AnnMonocots'),])
test_basalAnn<-subset(data_heatmap_percnbgene[ which(data_heatmap_percnbgene$hypothesis=='basalAnn'),])
test_Dicots<-subset(data_heatmap_percnbgene[ which(data_heatmap_percnbgene$hypothesis=='Dicots'),])

lm_AnnMonocots<-lm(test_AnnMonocots$scaffold_length~test_AnnMonocots$perc_nb_genes)
summary(lm_AnnMonocots)


p_AnnMonocots <- paste("p=", format(round(lmp(lm_AnnMonocots), 2), nsmall = 2), collapse = NULL)
lm_basalAnn<-lm(test_basalAnn$scaffold_length~test_basalAnn$perc_nb_genes)
summary(lm_basalAnn)
p_basalAnn <- paste("p=", format(round(lmp(lm_basalAnn), 2), nsmall = 2), collapse = NULL)
lm_Dicots<-lm(test_Dicots$scaffold_length~test_Dicots$perc_nb_genes)
summary(lm_Dicots)
p_Dicots <- paste("p=", format(round(lmp(lm_Dicots), 2), nsmall = 2), collapse = NULL)


pdf(file="scaffold_lentgh-vs-perc_nb_genes.pdf", width=12, height=12)
# inset
a <- ggplot(data_heatmap_percnbgene, aes(x=factor(hypothesis), y=perc_nb_genes, fill=hypothesis)) + geom_boxplot(outlier.size=1.5, notch=TRUE)+ theme(axis.title.x=element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), plot.background = element_rect(fill = "white"),legend.position = "none")
# large graph
#with trend lines# b <- ggplot(data_heatmap_percnbgene, aes(x=scaffold_length,, y=perc_nb_genes, colour=hypothesis))+ geom_point() + xlab("Scaffold length") + ylab("relative number of genes supporting a topology") + geom_rug(position="jitter", size=.2) + theme(legend.position = c(0.79,0.7), legend.title=element_blank(), legend.justification = c(-0.1,1)) + theme(legend.key=element_blank()) + stat_smooth(method=lm, se=FALSE)
#with legend# b <- ggplot(data_heatmap_percnbgene, aes(x=scaffold_length,, y=perc_nb_genes, colour=hypothesis))+ geom_point() + xlab("Scaffold length") + ylab("relative number of genes supporting a topology") + geom_rug(position="jitter", size=.2) + theme(legend.position = c(0.79,0.7), legend.title=element_blank(), legend.justification = c(-0.1,1)) + theme(legend.key=element_blank())
b <- ggplot(data_heatmap_percnbgene, aes(x=scaffold_length, y=perc_nb_genes, colour=hypothesis)) + geom_point() + geom_rug(position="jitter", size=.2) + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = "none", legend.title=element_blank(), legend.justification = c(-0.1,1)) + theme(legend.key=element_blank()) + stat_smooth(method=lm, se=FALSE) +
  annotate("text", x=20000, y=10, label= p_AnnMonocots, colour="red") +
  annotate("text", x=20000, y=38, label= p_basalAnn, colour="darkgreen") +
  annotate("text", x=20000, y=4, label= p_Dicots, colour="blue")
  

## left-sided plot
plot_left <- ggplot(data_heatmap_percnbgene, aes(perc_nb_genes, fill=hypothesis)) + geom_density(alpha=.5) + coord_flip() + scale_fill_manual(values = c("red", "green", "blue")) + theme(legend.position = "none", axis.text.y = element_blank()) + ylab("density") + xlab("relative number of genes supporting a topology (%)") + scale_y_reverse()
density_length <- ggplot(data_heatmap_percnbgene, aes(scaffold_length)) + geom_line(stat="density") + geom_density(fill="black", colour=NA, alpha=.2) + scale_fill_manual(values = c("red", "green", "blue")) + theme(legend.position = "none", axis.text.y = element_blank()) + ylab("density") + xlab("Scaffold length (kb)") 

blankPlot <- ggplot()+geom_blank(aes(1,1))+theme(plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank())


c <- grid.arrange(plot_left, b, blankPlot, density_length, nrow=2, ncol=2, widths=c(1,4), heights=c(4, 1))

#grid.newpage()
vpb_ <- viewport(width = 1, height = 1, x = 0.5, y = 0.5)  # the larger map
vpa_ <- viewport(width = 0.4, height = 0.3, x = 0.8, y = 0.847)  # the inset in upper right
# print(b, vp = vpb_)
print(c, vp = vpb_)
print(a, vp = vpa_)


dev.off()


#### making an histogram with frequency of data


class(data)<- "numeric"
plot(data$scaffold_length, data$nb_gene_background)

bins<-seq(0, 0.4, length.out = 100)
hist_percgene_background<-subset(data, select= perc_gene_background)
hist(hist_percgene_background, breaks = bins)
hist_percgene_Dicots<-subset(data, select= perc_gene_Dicots)
hist(hist_percgene_Dicots, breaks = bins)
hist_percgene_AnnMonocots<-subset(data, select= perc_gene_AnnMonocots)
hist(hist_percgene_background, breaks =bins)
hist_percgene_basalAnnona<-subset(data, select= perc_gene_basalAnnona)
hist(hist_percgene_basalAnnona, breaks = bins)
