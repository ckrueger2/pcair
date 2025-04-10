library(SNPRelate)
library(GENESIS)
library(GWASTools)
library(data.table)
library(dplyr)
library(readxl)
library(knitr)
library(ggplot2)
library(argparse)

#Rscript ~/Documents/quant_bio/03pcair.R --BfileDir ~/Documents/quant_bio/ -o ~/Documents/quant_bio/pcair/

#LD pruned plink files to GDS files
snpgdsBED2GDS(bed.fn = file.path("~/Documents/quant_bio/ld_pruned_all_hg38.bed"), 
              bim.fn = file.path("~/Documents/quant_bio/ld_pruned_all_hg38.bim"), 
              fam.fn = file.path("~/Documents/quant_bio/ld_pruned_all_hg38.fam"),
              cvt.chr = "char",
              out.gdsfn = file.path("~/Documents/quant_bio/pcair/merged_races.gds"))

#Close any open GDS files from previous runs
if (exists("gdsfile") && inherits(gdsfile, "SNPGDSFileClass")) {
  snpgdsClose(gdsfile)  # Properly close the previous GDS file
  rm(gdsfile)
  cat("Closed existing GDS file.\n")
}
if (exists("geno") && inherits(geno, "GdsGenotypeReader")) {
  close(geno)
  rm(geno)
  cat("Closed existing GenotypeData reader.\n")
}

#Define GDS file
gdsfile_path <- file.path("~/Documents/quant_bio/pcair/merged_races.gds")

#Open gdsfile
gdsfile <- snpgdsOpen(gdsfile_path)

#Create KING matrix
king_result <- snpgdsIBDKING(gdsfile)
king_matrix <- king_result$kinship
king_matrix[1:5,1:5] #Check old headers

#Rename matrix columns and check via printing first 5 columns
colnames(king_matrix) <- king_result$sample.id
row.names(king_matrix) <- king_result$sample.id
king_matrix[1:5,1:5] #Check new headers

#Close the SNPRelate 
snpgdsClose(gdsfile)

#Open the gdsfile as GdsGenotypeReader
geno <- GdsGenotypeReader(filename = file.path("~/Documents/quant_bio/pcair/merged_races.gds"))
genoData <- GenotypeData(geno)

#Run PC-AiR on pruned SNPs
pcair_run <- pcair(genoData, kinobj = king_matrix, divobj = king_matrix)

#Plot PCA-Air PCs
pdf(file.path("~/Documents/quant_bio/pcair/", "PC-AiR_plots.pdf"))

plot(pcair_run, vx = 1, vy = 2) # plot PCs 1 and 2
plot(pcair_run, vx = 1, vy = 3) # plot PCs 1 and 3
plot(pcair_run, vx = 1, vy = 4) # plot PCs 1 and 4
plot(pcair_run, vx = 2, vy = 3) # plot PCs 2 and 3
plot(pcair_run, vx = 3, vy = 4) # plot PCs 3 and 4

dev.off()

#Save eigenvalues
value_file <- file.path("~/Documents/quant_bio/pcair/pca.eigenval") #Save output path
eigenval <- pcair_run$values   # Eigenvalues for each PC
eigenval_subset <- eigenval[1:10, drop = FALSE] #Keep only first 10 PCs
write.table(eigenval_subset, file.path("~/Documents/quant_bio/pcair/pca.eigenval"), quote = FALSE, row.names = FALSE, col.names = FALSE)

#Save eigenvectors
#Get the summary of PC-AiR results
pcair_summary <- summary(pcair_run)

#Extract the eigenvectors (PC scores)
eigenvec <- pcair_summary$vectors

#Prepare the data frame for eigenvectors
fam_data <- read.table(file.path("~/Documents/quant_bio/", "ld_pruned_all_hg38.fam"), header = FALSE)
colnames(fam_data) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
eigenvec_df <- data.frame(FID = fam_data$FID, IID = fam_data$IID, eigenvec)

#Save eigenvectors to directory
eigenvec_subset <- eigenvec_df[1:12] #Keep only first 10 PCs
colnames(eigenvec_subset) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
eigenvec_path <- file.path("~/Documents/quant_bio/pcair/pca.eigenvec")  # Specify your output path
write.table(eigenvec_subset, eigenvec_path, quote = FALSE, row.names = FALSE, col.names = TRUE)

#Read in PCA files
pops = fread(file.path("~/Documents/quant_bio/all_hg38.psam"))
evec = fread(file.path("~/Documents/quant_bio/pcair/pca.eigenvec"))
eval = fread(file.path("~/Documents/quant_bio/pcair/pca.eigenval"))

all = left_join(evec, pops, by = c("IID" = "#IID"))

#Plot PCA
pdf(file.path("~/Documents/quant_bio/pcair/", "pca_plots.pdf"))

ggplot(all,aes(x=PC1,y=PC2,color=Population,label=SuperPop)) + geom_text(size=1.25)+ theme_bw(16)
ggplot(all,aes(x=PC1,y=PC3,color=Population,label=SuperPop)) + geom_text(size=1.25)+ theme_bw(16)
ggplot(all,aes(x=PC1,y=PC4,color=Population,label=SuperPop)) + geom_text(size=1.25) + theme_bw(16)
ggplot(all,aes(x=PC2,y=PC3,color=Population,label=SuperPop)) + geom_text(size=1.25)+ theme_bw(16)
ggplot(all,aes(x=PC3,y=PC4,color=Population,label=SuperPop)) + geom_text(size=1.25) + theme_bw(16)

dev.off()

#Use the dplyr function mutate to add a column named PC of numbers 1-10 and proportion variance explained by each PC
eval = mutate(eval,PC=as.factor(1:10),pve=V1/sum(V1))

#Plot PCs
pdf(file.path("~/Documents/quant_bio/pcair/", "scree_plot.pdf"))
ggplot(eval, aes(x=PC, y=pve, group=1)) + geom_point() + geom_line() + ylab("variance explained")+ theme_bw(16)
dev.off()