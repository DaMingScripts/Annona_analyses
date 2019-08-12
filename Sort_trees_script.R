library(PhySortR)

setwd("~/Dropbox/HUB/4._Annona_Genome/RESULTS/OrthoMCL")

### count number of trees supporting a dicot group
dicot_0=sortTrees(target.groups = "Amur,Atha", in.dir= "4_species_orthologs_longer100aa_aligned/", out.dir="support_eudicots/", mode = "c", min.support = 0)
dicot_50=sortTrees(target.groups = "Amur,Atha", in.dir= "4_species_orthologs_longer100aa_aligned/", out.dir="support_eudicots/BP50/", mode = "c", min.support = 0.50)
dicot_70=sortTrees(target.groups = "Amur,Atha", in.dir= "4_species_orthologs_longer100aa_aligned/", out.dir="support_eudicots/BP70/", mode = "c", min.support = 0.70)
dicot_95=sortTrees(target.groups = "Amur,Atha", in.dir= "4_species_orthologs_longer100aa_aligned/", out.dir="support_eudicots/BP95/", mode = "c", min.support = 0.95)


### count the number of trees supporting a Annona-monocot group
AnnMonocot_0=sortTrees(target.groups = "Amur,Osat", in.dir= "4_species_orthologs_longer100aa_aligned/", out.dir="support_AnnMonocots/", mode = "c", min.support = 0)
AnnMonocot_50=sortTrees(target.groups = "Amur,Osat", in.dir= "4_species_orthologs_longer100aa_aligned/", out.dir="support_AnnMonocots/BP50/", mode = "c", min.support = 0.50)
AnnMonocot_70=sortTrees(target.groups = "Amur,Osat", in.dir= "4_species_orthologs_longer100aa_aligned/", out.dir="support_AnnMonocots/BP70/", mode = "c", min.support = 0.70)
AnnMonocot_95=sortTrees(target.groups = "Amur,Osat", in.dir= "4_species_orthologs_longer100aa_aligned/", out.dir="support_AnnMonocots/BP95/", mode = "c", min.support = 0.95)


### count the number of trees supporting a basal Annona
basal_0=sortTrees(target.groups = "Osat,Atha", in.dir= "4_species_orthologs_longer100aa_aligned/", out.dir="support_basal_annona/", mode = "l", min.support = 0)
basal_50=sortTrees(target.groups = "Atha,Osat", in.dir= "4_species_orthologs_longer100aa_aligned/", out.dir="support_basal_annona/BP50/", mode = "l", min.support = 0.50)
basal_70=sortTrees(target.groups = "Atha,Osat", in.dir= "4_species_orthologs_longer100aa_aligned/", out.dir="support_basal_annona/BP70/", mode = "l", min.support = 0.70)
basal_95=sortTrees(target.groups = "Atha,Osat", in.dir= "4_species_orthologs_longer100aa_aligned/", out.dir="support_basal_annona/BP95/", mode = "l", min.support = 0.95)
