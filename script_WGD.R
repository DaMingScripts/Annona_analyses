# for the MacBook Air
# setwd("/Volumes/Dropbox/Dropbox/HUB/4._Annona_Genome/RESULTS/WGD_analysis")

# For the office desktop
setwd("~/Dropbox/HUB/4._Annona_Genome/RESULTS/WGD_analysis")
# For the Ubuntu laptop
setwd("/media/hin175/Professional/Dropbox/HUB/4._Annona_Genome/RESULTS/WGD_analysis")


library(mixtools)

data<-read.csv("KS_annona.csv")
data2<-as.matrix(data)
#hist_data<-hist(data2, breaks=1000,xlim=c(0,8.0), xlab = "Ks", main = "Age distribution of the A. muricata paralogs based on K S values")

data_Annona_Annona <- subset(data2, select = Annona_Annona)
data_Nelumbo_Arabidopsis <- subset(data2, select = Nnu_Ath)
data_Amborella_Oryza <- subset(data2, select = Atr_Osat)
data_Annona_Oryza <- subset(data2, select = Asq_Osa)
data_Annona_Amborella <- subset(data2, select = Asq_Atr)
data_Annona_Arabidopsis <- subset(data2, select = Asq_Ath)
data_Annona_Nelumbo <- subset(data2, select = Asq_Nnu)

#pdf(file="Ks_Annona_genome.pdf", width="800", height="300")

op <- par(mfrow = c(4, 2), pty = "s", ann=TRUE, mar=c(0.5, 0, 1.5,0), oma=c(2,0,0,0))

#layout(matrix(1:6, 3, 2))
#layout.show(6)
hist_data_Annona_Nelumbo<-hist(data_Annona_Nelumbo, breaks=1000,xlim=c(0,8.0), xlab = "", axes=FALSE, main = "Annona/Nelumbo")
axis(side=2)
hist_data_Nelumbo_Arabidopsis<-hist(data_Nelumbo_Arabidopsis, breaks=1000,xlim=c(0,8.0), axes=FALSE, xlab = "", main = "Nelumbo/Arabidopsis")
axis(side=2)
hist_data_Amborella_Oryza<-hist(data_Amborella_Oryza, breaks=1000,xlim=c(0,8.0), xlab = "", axes=FALSE, main = "Amborella/Oryza")
axis(side=2)
hist_data_Annona_Oryza<-hist(data_Annona_Oryza, breaks=1000,xlim=c(0,8.0), xlab = "", axes=FALSE, main = "Annona/Oryza")
axis(side=2)
hist_data_Annona_Amborella<-hist(data_Annona_Amborella, breaks=1000,xlim=c(0,8.0), xlab = "", axes=FALSE, main = "Annona/Amborella")
axis(side=2)
hist_data_Annona_Arabidopsis<-hist(data_Annona_Arabidopsis, breaks=1000,xlim=c(0,8.0), xlab = "", axes=FALSE, main = "Annona/Arabidopsis")
axis(side=2)
axis(side=1)
hist_data_Annona_Annona<-hist(data_Annona_Annona, breaks=1000,xlim=c(0,8.0), axes=TRUE, xlab="", main = "A. muricata paralogs")
par(op)

par(op)
#dev.off()


# wait = faithful$waiting
data_mix_AnnAnn = as.numeric(hist_data_Annona_Annona$counts)
data_mix_AnnNel = as.numeric(hist_data_Annona_Nelumbo$counts)
data_mix_NelAra = as.numeric(hist_data_Nelumbo_Arabidopsis$counts)
data_mix_AmbOry = as.numeric(hist_data_Amborella_Oryza$counts)
data_mix_AnnOry = as.numeric(hist_data_Annona_Oryza$counts)
data_mix_AnnAmb = as.numeric(hist_data_Annona_Amborella$counts)
data_mix_AnnAnn = as.numeric(hist_data_Annona_Arabidopsis$counts)
#K = 2
#start.par <- mean(wait2, na.rm = TRUE) + sd(wait2, na.rm = TRUE) * runif(K)

data_mix_filtered <- data2[1:63]

mixmdl = normalmixEM(data_mix_filtered, k=3)
plot(mixmdl,which=2, xlim=c(0,8.0))
#plot(hist_data, breaks=1000,xlim=c(0,8.0), xlab = "Ks", main = "Age distribution of the A. muricata paralogs based on K S values")
lines(density(data_mix_AnnAnn), lty=2, lwd=2)
