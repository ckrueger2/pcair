library(dplyr)
library(data.table)

#PC-AIR
#Read in PCA files
pops_pcair = fread(file.path("~/Documents/quant_bio/all_hg38.psam"))
evec_pcair = fread(file.path("~/Documents/quant_bio/pcair/pcair.eigenvec"))
eval_pcair = fread(file.path("~/Documents/quant_bio/pcair/pcair.eigenval"))

all_pcair = left_join(evec_pcair, pops_pcair, by = c("IID" = "#IID"))

#Plot PCA
pdf(file.path("~/Documents/quant_bio/pcair/", "pcair_pca_plots.pdf"))

ggplot(all_pcair,aes(x=PC1,y=PC2,color=Population,label=SuperPop)) + geom_text(size=1.25)+ theme_bw(16)
ggplot(all_pcair,aes(x=PC2,y=PC3,color=Population,label=SuperPop)) + geom_text(size=1.25)+ theme_bw(16)
ggplot(all_pcair,aes(x=PC3,y=PC4,color=Population,label=SuperPop)) + geom_text(size=1.25) + theme_bw(16)
ggplot(all_pcair,aes(x=PC4,y=PC5,color=Population,label=SuperPop)) + geom_text(size=1.25) + theme_bw(16)

dev.off()

#Use the dplyr function mutate to add a column named PC of numbers 1-10 and proportion variance explained by each PC
eval_pcair = mutate(eval_pcair,PC=as.factor(1:10),pve=V1/sum(V1))

#Plot individual scree plot
pdf(file.path("~/Documents/quant_bio/pcair/", "pcair_scree_plot.pdf"))
ggplot(eval_pcair, aes(x=PC, y=pve, group=1)) + geom_point() + geom_line() + ylab("variance explained")+ theme_bw(16)
dev.off()

#STANDARD_PCA
#Read in PCA files
pops_pca = fread(file.path("~/Documents/quant_bio/all_hg38.psam"))
evec_pca = fread(file.path("~/Documents/quant_bio/pcair/pca.eigenvec"))
eval_pca = fread(file.path("~/Documents/quant_bio/pcair/pca.eigenval"))

all_pca = left_join(evec_pca, pops_pca, by = c("IID" = "#IID"))

#Plot PCA
pdf(file.path("~/Documents/quant_bio/pcair/", "pca_standard_pca_plots.pdf"))

ggplot(all_pca,aes(x=PC1,y=PC2,color=Population,label=SuperPop)) + geom_text(size=1.25)+ theme_bw(16)
ggplot(all_pca,aes(x=PC2,y=PC3,color=Population,label=SuperPop)) + geom_text(size=1.25)+ theme_bw(16)
ggplot(all_pca,aes(x=PC3,y=PC4,color=Population,label=SuperPop)) + geom_text(size=1.25) + theme_bw(16)
ggplot(all_pca,aes(x=PC4,y=PC5,color=Population,label=SuperPop)) + geom_text(size=1.25) + theme_bw(16)

dev.off()

#Use the dplyr function mutate to add a column named PC of numbers 1-10 and proportion variance explained by each PC
eval_pca = mutate(eval_pca,PC=as.factor(1:10),pve=V1/sum(V1))

#Plot individual scree plot
pdf(file.path("~/Documents/quant_bio/pcair/", "pca_standard_scree_plot.pdf"))
ggplot(eval_pca, aes(x=PC, y=pve, group=1)) + geom_point() + geom_line() + ylab("variance explained")+ theme_bw(16)
dev.off()