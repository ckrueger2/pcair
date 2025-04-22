library(SNPRelate)
library(GENESIS)
library(GWASTools)
library(data.table)
library(dplyr)
library(readxl)
library(knitr)
library(ggplot2)
library(argparse)
library(vegan)

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

#PC-AIR ANALYSIS
#Create KING matrix
king_result <- snpgdsIBDKING(gdsfile)
king_matrix <- king_result$kinship
king_matrix[1:5,1:5] #Check old headers

#Rename matrix columns and check via printing first 5 columns
colnames(king_matrix) <- king_result$sample.id
row.names(king_matrix) <- king_result$sample.id
king_matrix[1:5,1:5] #Check new headers

#Save king results
save(king_result, king_matrix, file = paste0(results_base, "_king.RData"))

#Close the SNPRelate 
snpgdsClose(gdsfile)

#Open the gdsfile as GdsGenotypeReader
geno <- GdsGenotypeReader(filename = file.path("~/Documents/quant_bio/pcair/merged_races.gds"))
genoData <- GenotypeData(geno)

#Run PC-AiR on pruned SNPs
pcair_run <- pcair(genoData, kinobj = king_matrix, divobj = king_matrix)

#Close the GenotypeData reader
close(geno)

#Save PC-AIR results
save(pcair_run, file = pcair_results)

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
eigenval <- pcair_run$values   #Eigenvalues for each PC
eigenval_subset <- eigenval[1:10, drop = FALSE] #Keep only first 10 PCs
write.table(eigenval_subset, file.path("~/Documents/quant_bio/pcair/pcair.eigenval"), quote = FALSE, row.names = FALSE, col.names = FALSE)

#Save eigenvectors
#Get the summary of PC-AiR results
pcair_summary <- summary(pcair_run)

#Extract the eigenvectors (PC scores)
eigenvec <- pcair_summary$vectors

#Prepare the data frame for eigenvectors
fam_data <- read.table(file.path("~/Documents/quant_bio/", "ld_pruned_all_hg38.fam"), header = FALSE)
colnames(fam_data) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
eigenvec_df_pcair <- data.frame(FID = fam_data$FID, IID = fam_data$IID, eigenvec)

#Save eigenvectors to directory
eigenvec_subset_pcair <- eigenvec_df_pcair[1:12] #Keep only first 10 PCs
colnames(eigenvec_subset_pcair) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
eigenvec_path_pcair <- file.path("~/Documents/quant_bio/pcair/pcair.eigenvec")  # Specify your output path
write.table(eigenvec_subset_pcair, eigenvec_path_pcair, quote = FALSE, row.names = FALSE, col.names = TRUE)

#Save results
save(pcair_run, fam_data, pcair_summary, eigenvec_df_pcair, file = paste0(results_base, "_pcair_vectors.RData"))

#STANDARD PCA
#Re-open GDS file
gdsfile <- snpgdsOpen(gdsfile_path)

#Run standard PCA
pca_run <- snpgdsPCA(gdsfile, num.thread = 2)

#Save PCA results
save(pca_run, file = pca_results)

#Save eigenvalues for standard PCA
eigenval_pca <- pca_run$eigenval
eigenval_subset_pca <- eigenval_pca[1:10, drop = FALSE]
write.table(eigenval_subset_pca, file.path("~/Documents/quant_bio/pcair/pca.eigenval"), quote = FALSE, row.names = FALSE, col.names = FALSE)

#Save eigenvectors for standard PCA
eigenvec_pca <- pca_run$eigenvect[, 1:10]
sample_ids <- read.gdsn(index.gdsn(gdsfile, "sample.id"))
eigenvec_df_pca <- data.frame(FID = fam_data$FID, IID = fam_data$IID, eigenvec_pca)
colnames(eigenvec_df_pca) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
write.table(eigenvec_df_pca, file.path("~/Documents/quant_bio/pcair/pca.eigenvec"), quote = FALSE, row.names = FALSE, col.names = TRUE)

#Close GDS file
snpgdsClose(gdsfile)

#Save results
save(pca_run, eigenvec_df_pca, file = paste0(results_base, "_pca_vectors.RData"))