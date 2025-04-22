library(dplyr)
library(data.table)

#COMBINE DATA
#Read in PCA and PC-AIR files
pops = fread(file.path("~/Documents/quant_bio/all_hg38.psam"))
evec_pca = fread(file.path("~/Documents/quant_bio/pcair/pca.eigenvec"))
evec_pcair = fread(file.path("~/Documents/quant_bio/pcair/pcair.eigenvec"))
eval_pca = fread(file.path("~/Documents/quant_bio/pcair/pca.eigenval"))
eval_pcair = fread(file.path("~/Documents/quant_bio/pcair/pcair.eigenval"))

#Join with population info
all_pca = left_join(evec_pca, pops, by = c("IID" = "#IID"))
all_pcair = left_join(evec_pcair, pops, by = c("IID" = "#IID"))

#Add a method column
all_pca$Method = "Standard PCA"
all_pcair$Method = "PC-AIR"

#Combine datasets
all_combined = rbind(all_pca, all_pcair)

#Save combined results
save(pops, evec_pca, evec_pcair, eval_pca, eval_pcair, all_pca, all_pcair, all_combined, file = combined_results)

#Plot comparison of PC1 vs PC2 for both methods
pdf(file.path("~/Documents/quant_bio/pcair/", "pca_comparison_plots.pdf"))

#PC1 vs PC2 by Population
ggplot(all_combined, aes(x=PC1, y=PC2, color=Population, shape=Method)) +
  geom_point(alpha=0.7) +
  theme_bw(16) +
  facet_wrap(~Method, scales="free") +
  labs(title="Comparison of Standard PCA and PC-AIR")

#PC1 vs PC2 by SuperPop
ggplot(all_combined, aes(x=PC1, y=PC2, color=SuperPop, shape=Method)) +
  geom_point(alpha=0.7) +
  theme_bw(16) +
  facet_wrap(~Method, scales="free") +
  labs(title="Comparison of Standard PCA and PC-AIR by Super Population")

dev.off()

#PERMANOVA
#Extract PC coordinates for both methods
pca_coords <- all_pca[, c("PC1", "PC2", "PC3", "PC4", "PC5")]
pcair_coords <- all_pcair[, c("PC1", "PC2", "PC3", "PC4", "PC5")]

#Combined dataset of method and population
permanova_data <- rbind(
  data.frame(
    Method = "Standard PCA",
    Population = all_pca$Population,
    SuperPop = all_pca$SuperPop,
    PC1 = pca_coords$PC1,
    PC2 = pca_coords$PC2,
    PC3 = pca_coords$PC3,
    PC4 = pca_coords$PC4,
    PC5 = pca_coords$PC5
  ),
  data.frame(
    Method = "PC-AIR",
    Population = all_pcair$Population,
    SuperPop = all_pcair$SuperPop,
    PC1 = pcair_coords$PC1,
    PC2 = pcair_coords$PC2,
    PC3 = pcair_coords$PC3,
    PC4 = pcair_coords$PC4,
    PC5 = pcair_coords$PC5
  )
)

#Matrix of PC data for PERMANOVA
pc_matrix <- as.matrix(permanova_data[, c("PC1", "PC2", "PC3", "PC4", "PC5")])

#Run PERMANOVA to test for differences between methods
perm_method <- adonis2(
  formula = pc_matrix ~ Method, 
  data = permanova_data, 
  method = "euclidean",
  permutations = 999
)

#Save results
save(pca_coords, pcair_coords, permanova_data, pc_matrix, perm_method, file = paste0(results_base, "_permanova_method.RData"))

#Run PERMANOVA for differences between populations and method combo
perm_pop_method <- adonis2(
  formula = pc_matrix ~ Population + Method, 
  data = permanova_data, 
  method = "euclidean",
  permutations = 999
)

#Save results
save(perm_pop_method, file = paste0(results_base, "_permanova_popmethod.RData"))

#Run PERMANOVA to test for interaction between method and population
perm_interaction <- adonis2(
  formula = pc_matrix ~ Population * Method, 
  data = permanova_data, 
  method = "euclidean",
  permutations = 999
)

#Save results
save(perm_interaction, file = paste0(results_base, "_permanova_interaction.RData"))

#Run PERMANOVA by super population
perm_superpop <- adonis2(
  formula = pc_matrix ~ SuperPop * Method, 
  data = permanova_data, 
  method = "euclidean",
  permutations = 999
)

#Save PERMANOVA results
save(pca_coords, pcair_coords, permanova_data, perm_method, perm_pop_method, perm_interaction, perm_superpop, file = permanova_results)

#COMPARE SCREE PLOTS
eval_pca = mutate(eval_pca, PC=as.factor(1:10), pve=V1/sum(V1), Method="Standard PCA")
eval_pcair = mutate(eval_pcair, PC=as.factor(1:10), pve=V1/sum(V1), Method="PC-AIR")
eval_combined = rbind(eval_pca, eval_pcair)

#Save scree plot data
save(eval_pca, eval_pcair, eval_combined, file = scree_results)

#Plot
pdf(file.path("~/Documents/quant_bio/pcair/", "scree_plot_comparison.pdf"))
ggplot(eval_combined, aes(x=PC, y=pve, group=Method, color=Method)) + 
  geom_point() + 
  geom_line() + 
  ylab("Proportion of Variance Explained") + 
  theme_bw(16) +
  labs(title="Comparison of Variance Explained: Standard PCA vs PC-AIR")
dev.off()

#COMPARE DISTANCE MATRICIES WITH MANTEL TEST
#Calculate distance matrices for both methods
pca_dist <- dist(pca_coords, method="euclidean")
pcair_dist <- dist(pcair_coords, method="euclidean")

#Save distance matrices
save(pca_dist, pcair_dist, file = paste0(results_base, "_distance_matrices.RData"))

#Compare the correlation between distance matrices
dist_correlation <- mantel(pca_dist, pcair_dist, method="spearman", permutations=999)

#Save results
save(dist_correlation, file = paste0(results_base, "_mantel.RData"))

#Save the Mantel test
sink(file.path("~/Documents/quant_bio/pcair/", "mantel_test_result.txt"))
cat("Mantel test comparing distance matrices from Standard PCA and PC-AIR:\n")
print(dist_correlation)
sink()

#NMDS plots
pca_nmds <- metaMDS(pca_dist, k=2)
save(pca_nmds, file = paste0(results_base, "_pca_nmds.RData"))

pcair_nmds <- metaMDS(pcair_dist, k=2)
save(pca_nmds, pcair_nmds, file = paste0(results_base, "_nmds.RData"))

#Data frames for NMDS plots
pca_nmds_df <- data.frame(
  NMDS1 = pca_nmds$points[,1],
  NMDS2 = pca_nmds$points[,2],
  Population = all_pca$Population,
  SuperPop = all_pca$SuperPop,
  Method = "Standard PCA"
)

pcair_nmds_df <- data.frame(
  NMDS1 = pcair_nmds$points[,1],
  NMDS2 = pcair_nmds$points[,2],
  Population = all_pcair$Population,
  SuperPop = all_pcair$SuperPop,
  Method = "PC-AIR"
)

nmds_combined <- rbind(pca_nmds_df, pcair_nmds_df)

save(pca_dist, pcair_dist, dist_correlation, pca_nmds, pcair_nmds, pca_nmds_df, pcair_nmds_df, nmds_combined, file = distance_results)

#Plot NMDS 
pdf(file.path("~/Documents/quant_bio/pcair/", "nmds_comparison.pdf"))
ggplot(nmds_combined, aes(x=NMDS1, y=NMDS2, color=Population, shape=Method)) + 
  geom_point(alpha=0.7) + 
  theme_bw(16) +
  facet_wrap(~Method) +
  labs(title="NMDS Comparison of Standard PCA and PC-AIR")

ggplot(nmds_combined, aes(x=NMDS1, y=NMDS2, color=SuperPop, shape=Method)) + 
  geom_point(alpha=0.7) + 
  theme_bw(16) +
  facet_wrap(~Method) +
  labs(title="NMDS Comparison of Standard PCA and PC-AIR by Super Population")
dev.off()

#FINAL SUMMARY
summary_stats <- permanova_data %>%
  group_by(Method, Population) %>%
  summarize(
    PC1_mean = mean(PC1),
    PC1_sd = sd(PC1),
    PC2_mean = mean(PC2),
    PC2_sd = sd(PC2)
  )

#Save summary statistics
save(summary_stats, file = summary_results)
write.csv(summary_stats, file.path("~/Documents/quant_bio/pcair/", "pc_summary_by_population.csv"), row.names = FALSE)

#Comparison
pdf(file.path("~/Documents/quant_bio/pcair/", "method_comparison_report.pdf"))

#Plot density distributions of PC1 and PC2 by method
ggplot(permanova_data, aes(x=PC1, fill=Method)) + 
  geom_density(alpha=0.5) + 
  theme_bw(16) +
  labs(title="Distribution of PC1 by Method")

ggplot(permanova_data, aes(x=PC2, fill=Method)) + 
  geom_density(alpha=0.5) + 
  theme_bw(16) +
  labs(title="Distribution of PC2 by Method")

#Plot PC1 vs PC2 by method and super population (all in one plot)
ggplot(permanova_data, aes(x=PC1, y=PC2, color=SuperPop, shape=Method)) + 
  geom_point(alpha=0.7) + 
  theme_bw(16) +
  labs(title="PC1 vs PC2 by Method and Super Population")

dev.off()

#Clean up files
if (exists("gdsfile") && inherits(gdsfile, "SNPGDSFileClass")) {
  snpgdsClose(gdsfile)
  rm(gdsfile)
  cat("Closed existing GDS file.\n")
}
if (exists("geno") && inherits(geno, "GdsGenotypeReader")) {
  close(geno)
  rm(geno)
  cat("Closed existing GenotypeData reader.\n")
}