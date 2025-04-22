library(SNPRelate)
library(GENESIS)
library(GWASTools)
library(data.table)
library(dplyr)
library(ggplot2)
library(vegan)

#Define base directories
base_dir <- "~/Documents/quant_bio/"
pcair_dir <- file.path(base_dir, "pcair_by_population")

#Create main directory
dir.create(pcair_dir, recursive = TRUE)

#Read population information
pops <- fread(file.path(base_dir, "all_hg38.psam"))
fam_data <- read.table(file.path(base_dir, "ld_pruned_all_hg38.fam"), header = FALSE)
colnames(fam_data) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")

#Get list of all populations
populations <- unique(pops$Population)

#Define GDS file exists
gds_path <- file.path(base_dir, "pcair/merged_races.gds")
snpgdsBED2GDS(bed.fn = file.path(base_dir, "ld_pruned_all_hg38.bed"),
              bim.fn = file.path(base_dir, "ld_pruned_all_hg38.bim"),
              fam.fn = file.path(base_dir, "ld_pruned_all_hg38.fam"),
              cvt.chr = "char",
              out.gdsfn = gds_path)

#Define function to analyze by pop
run_population_analysis <- function(population) {
  
  #Create directory for this population
  pop_dir <- file.path(pcair_dir, population)
  dir.create(pop_dir, recursive = TRUE)
  
  #Define base results file for this population
  results_base <- file.path(pop_dir, paste0(population, "_analysis_results"))
  
  #Filter data for this population
  pop_samples <- pops[Population == population, `#IID`]
  
  #PC-AIR
  pcair_results <- paste0(results_base, "_pcair.RData")

  #Open GDS file
  gdsfile <- snpgdsOpen(gds_path)
  
  #Get all sample IDs
  all_samples <- read.gdsn(index.gdsn(gdsfile, "sample.id"))
  
  #Find which samples from this population exist in the GDS file
  valid_samples <- intersect(all_samples, pop_samples)
  
  if (length(valid_samples) == 0) {
    gds_samples_clean <- gsub("^[A-Za-z]+", "", all_samples)
    pop_samples_clean <- gsub("^[A-Za-z]+", "", pop_samples)
    
    matches <- match(pop_samples_clean, gds_samples_clean)
    valid_indices <- which(!is.na(matches))
    
    if (length(valid_indices) > 0) {
      valid_samples <- all_samples[matches[valid_indices]]
    } else {
      valid_samples <- all_samples
    }
  }
  
  #KING matrix
  king_result <- snpgdsIBDKING(gdsfile, sample.id = valid_samples)
  king_matrix <- king_result$kinship
  colnames(king_matrix) <- king_result$sample.id
  row.names(king_matrix) <- king_result$sample.id
  
  #Close the SNPRelate file
  snpgdsClose(gdsfile)
  
  #Open as GdsGenotypeReader
  geno <- GdsGenotypeReader(filename = gds_path)
  genoData <- GenotypeData(geno)
  
  #Run PC-AiR
  pcair_run <- pcair(genoData, 
                     kinobj = king_matrix, 
                     divobj = king_matrix,
                     sample.include = valid_samples)
  
  #Close the GenotypeData reader
  close(geno)
  
  #Save PC-AIR results
  save(pcair_run, file = pcair_results)
  
  #Plot PC pairs
  pc_combinations <- list(c(1, 2), c(1, 3), c(2, 3))
  for (pc_pair in pc_combinations) {
    png(file.path(pop_dir, paste0(population, "_PC-AiR_PC", pc_pair[1], "_PC", pc_pair[2], ".png")), 
        width = 800, height = 600)
    plot(pcair_run, vx = pc_pair[1], vy = pc_pair[2])
    dev.off()
  }
  
  #Save eigenvalues
  eigenval <- pcair_run$values[1:10]
  write.table(eigenval, file.path(pop_dir, paste0(population, "_pcair.eigenval")), quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  #Save eigenvectors
  pcair_summary <- summary(pcair_run)
  eigenvec <- pcair_summary$vectors
  
  #Match with fam data
  sample_ids <- pcair_run$sample.id
  fam_subset <- fam_data[match(sample_ids, fam_data$IID), ]
  
  #Create eigenvector dataframe
  eigenvec_df_pcair <- data.frame(FID = fam_subset$FID, IID = sample_ids, eigenvec)
  colnames(eigenvec_df_pcair) <- c("FID", "IID", paste0("PC", 1:ncol(eigenvec)))
  eigenvec_subset_pcair <- eigenvec_df_pcair[, 1:(2+10)] # Keep only first 10 PCs + FID/IID
  write.table(eigenvec_subset_pcair, file.path(pop_dir, paste0(population, "_pcair.eigenvec")), quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  #STANDARD PCA
  pca_results <- paste0(results_base, "_pca.RData")
  load(pca_results)
  
  #Open the main GDS file
  gdsfile <- snpgdsOpen(gds_path)
  
  #Get all sample IDs
  all_samples <- read.gdsn(index.gdsn(gdsfile, "sample.id"))
  
  #Find which samples from this population exist in the GDS file
  valid_samples <- intersect(all_samples, pop_samples)
  
  if (length(valid_samples) == 0) {
    gds_samples_clean <- gsub("^[A-Za-z]+", "", all_samples)
    pop_samples_clean <- gsub("^[A-Za-z]+", "", pop_samples)

    matches <- match(pop_samples_clean, gds_samples_clean)
    valid_indices <- which(!is.na(matches))
    
    if (length(valid_indices) > 0) {
      valid_samples <- all_samples[matches[valid_indices]]
    } else {
      valid_samples <- all_samples
    }
  }
  
  #Run standard PCA
  pca_run <- snpgdsPCA(gdsfile, sample.id = valid_samples, num.thread = 2)
  
  #Save
  save(pca_run, file = pca_results)
  
  #Save eigenvalues
  eigenval_pca <- pca_run$eigenval[1:10]
  write.table(eigenval_pca, file.path(pop_dir, paste0(population, "_pca.eigenval")), quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  #Save eigenvectors
  eigenvec_pca <- pca_run$eigenvect[, 1:10]
  
  #Match sample IDs 
  sample_ids <- pca_run$sample.id
  fam_subset <- fam_data[match(sample_ids, fam_data$IID), ]
  
  #Create eigenvector dataframe
  eigenvec_df_pca <- data.frame(FID = fam_subset$FID, IID = sample_ids, eigenvec_pca)
  colnames(eigenvec_df_pca) <- c("FID", "IID", paste0("PC", 1:10))
  write.table(eigenvec_df_pca, file.path(pop_dir, paste0(population, "_pca.eigenvec")), quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  #Close GDS file
  snpgdsClose(gdsfile)
  
  #COMBINE AND PLOT
  combined_results <- paste0(results_base, "_combined.RData")
  load(combined_results)
  
  #Read in PCA and PC-AIR files
  evec_pca = fread(file.path(pop_dir, paste0(population, "_pca.eigenvec")))
  evec_pcair = fread(file.path(pop_dir, paste0(population, "_pcair.eigenvec")))
  eval_pca = fread(file.path(pop_dir, paste0(population, "_pca.eigenval")))
  eval_pcair = fread(file.path(pop_dir, paste0(population, "_pcair.eigenval")))
  
  #Filter population data
  pop_info <- pops[pops$Population == population, ]
  
  #Join with population information
  all_pca = left_join(evec_pca, pop_info, by = c("IID" = "#IID"))
  all_pcair = left_join(evec_pcair, pop_info, by = c("IID" = "#IID"))
  
  #Add a method column
  all_pca$Method = "Standard PCA"
  all_pcair$Method = "PC-AIR"
  
  #Combine
  all_combined = rbind(all_pca, all_pcair)
  
  #Save
  save(all_pca, all_pcair, all_combined, file = combined_results)
  
  #Plot comparison of PC1 vs PC2 for both methods
  p1 <- ggplot(all_combined, aes(x=PC1, y=PC2, color=SuperPop, shape=Method)) +
    geom_point(alpha=0.7) +
    theme_bw(12) +
    facet_wrap(~Method, scales="free") +
    ggtitle(paste0("Comparison for ", population))
  
  ggsave(file.path(pop_dir, paste0(population, "_pca_comparison.png")), p1, width = 10, height = 8)
  
  #PERMANOVA
  permanova_results <- paste0(results_base, "_permanova.RData")
  load(permanova_results)
  
  #Extract PCs
  pca_coords <- all_pca[, c("PC1", "PC2", "PC3", "PC4", "PC5")]
  pcair_coords <- all_pcair[, c("PC1", "PC2", "PC3", "PC4", "PC5")]
  
  #Combine dataset
  permanova_data <- rbind(
    data.frame(
      Method = "Standard PCA",
      SuperPop = all_pca$SuperPop,
      PC1 = pca_coords$PC1,
      PC2 = pca_coords$PC2,
      PC3 = pca_coords$PC3,
      PC4 = pca_coords$PC4,
      PC5 = pca_coords$PC5
    ),
    data.frame(
      Method = "PC-AIR",
      SuperPop = all_pcair$SuperPop,
      PC1 = pcair_coords$PC1,
      PC2 = pcair_coords$PC2,
      PC3 = pcair_coords$PC3,
      PC4 = pcair_coords$PC4,
      PC5 = pcair_coords$PC5
    )
  )
  
  #Matrix for PERMANOVA
  pc_matrix <- as.matrix(permanova_data[, c("PC1", "PC2", "PC3", "PC4", "PC5")])
  
  #Check if SuperPop has more than one level
  superpop_levels <- length(unique(permanova_data$SuperPop))
  
  #Run PERMANOVA tests with interactions - only if we have sufficient factor levels
  #Run Method-only PERMANOVA
  perm_method <- adonis2(pc_matrix ~ Method, data = permanova_data, 
                         method = "euclidean", permutations = 999)
  
  #Only run SuperPop models if we have multiple SuperPop levels
  if (superpop_levels > 1) {
    perm_superpop_method <- adonis2(pc_matrix ~ SuperPop + Method, data = permanova_data, 
                                    method = "euclidean", permutations = 999)
    
    #Interaction test
    perm_interaction <- adonis2(pc_matrix ~ SuperPop * Method, data = permanova_data, 
                                method = "euclidean", permutations = 999)
  } else {
    perm_superpop_method <- NULL
    perm_interaction <- NULL
  }
  
  #Save
  save(permanova_data, pc_matrix, perm_method, perm_superpop_method, perm_interaction, file = permanova_results)
  sink(file.path(pop_dir, paste0(population, "_permanova_results.txt")))

  
  if (superpop_levels > 1) {
    cat("\n2. Effect of SuperPop and Method:\n")
    print(perm_superpop_method)
    cat("\n3. Interaction between SuperPop and Method:\n")
    print(perm_interaction)
  } else {
    cat("\nNote: Only one SuperPop level found for this population, so SuperPop PERMANOVA tests were skipped.\n")
  }
  sink()
  
  #DISTANCE MATRIX ANALYSIS MANTEL TEST
  distance_results <- paste0(results_base, "_distance.RData")
  load(distance_results)
  
  #If pc_matrix doesn't exist from PERMANOVA section, recreate it
  if (!exists("pc_matrix") || !exists("permanova_data")) {
    pca_coords <- all_pca[, c("PC1", "PC2", "PC3", "PC4", "PC5")]
    pcair_coords <- all_pcair[, c("PC1", "PC2", "PC3", "PC4", "PC5")]
    
    permanova_data <- rbind(
      data.frame(
        Method = "Standard PCA",
        SuperPop = all_pca$SuperPop,
        as.data.frame(pca_coords)
      ),
      data.frame(
        Method = "PC-AIR",
        SuperPop = all_pcair$SuperPop,
        as.data.frame(pcair_coords)
      )
    )
    
    pc_matrix <- as.matrix(permanova_data[, c("PC1", "PC2", "PC3", "PC4", "PC5")])
  }
  
  #Split the data by method
  pca_indices <- which(permanova_data$Method == "Standard PCA")
  pcair_indices <- which(permanova_data$Method == "PC-AIR")
  
  #Calculate distance matrices for both methods
  pca_dist <- dist(pc_matrix[pca_indices,], method="euclidean")
  pcair_dist <- dist(pc_matrix[pcair_indices,], method="euclidean")
  
  #Compare the correlation between distance matrices using Mantel test
  if (length(pca_indices) == length(pcair_indices)) {
    dist_correlation <- mantel(pca_dist, pcair_dist, method="spearman", permutations=999)
    
    # Save the Mantel test result to text file
    sink(file.path(pop_dir, paste0(population, "_mantel_test_result.txt")))
    print(dist_correlation)
    sink()
  } else {
    dist_correlation <- NULL
  }
  
  # Save distance matrices and correlation
  save(pca_dist, pcair_dist, dist_correlation, file = distance_results)
  
  #NMDS ANALYSIS
  nmds_results <- paste0(results_base, "_nmds.RData")
  load(nmds_results)
  
  #Load distance matrices if not already in environment
  if (!exists("pca_dist") || !exists("pcair_dist")) {
    load(distance_results)
  }
  
  #Create NMDS plots based on the distance matrices
  pca_nmds <- metaMDS(pca_dist, k=2, trymax=100)
  pcair_nmds <- metaMDS(pcair_dist, k=2, trymax=100)
  
  #Create data frames
  pca_superpop <- all_pca$SuperPop[match(names(pca_dist), all_pca$IID)]
  if (all(is.na(pca_superpop))) {
    pca_superpop <- all_pca$SuperPop[1:nrow(as.matrix(pca_dist))]
  }
  
  pca_nmds_df <- data.frame(
    NMDS1 = pca_nmds$points[,1],
    NMDS2 = pca_nmds$points[,2],
    SuperPop = pca_superpop,
    Method = "Standard PCA"
  )
  
  pcair_superpop <- all_pcair$SuperPop[match(names(pcair_dist), all_pcair$IID)]
  if (all(is.na(pcair_superpop))) {
    # Try matching by row index if IID matching fails
    pcair_superpop <- all_pcair$SuperPop[1:nrow(as.matrix(pcair_dist))]
  }
  
  pcair_nmds_df <- data.frame(
    NMDS1 = pcair_nmds$points[,1],
    NMDS2 = pcair_nmds$points[,2],
    SuperPop = pcair_superpop,
    Method = "PC-AIR"
  )
  
  nmds_combined <- rbind(pca_nmds_df, pcair_nmds_df)
  
  #Save
  save(pca_nmds, pcair_nmds, pca_nmds_df, pcair_nmds_df, nmds_combined, file = nmds_results)
  
  #Plot
  p_nmds <- ggplot(nmds_combined, aes(x=NMDS1, y=NMDS2, color=SuperPop, shape=Method)) + 
    geom_point(alpha=0.7) + 
    theme_bw(12) +
    facet_wrap(~Method) +
    ggtitle(paste0("NMDS Comparison for ", population))
  
  ggsave(file.path(pop_dir, paste0(population, "_nmds_comparison.png")), 
         p_nmds, width = 10, height = 8)
  
  #SCREE PLOTS
  scree_results <- paste0(results_base, "_scree.RData")
  load(scree_results)
  
  #Load eigenvalues
  eval_pca = fread(file.path(pop_dir, paste0(population, "_pca.eigenval")))
  eval_pcair = fread(file.path(pop_dir, paste0(population, "_pcair.eigenval")))
  
  #Calculate proportion of variance explained
  eval_pca = mutate(eval_pca, PC=as.factor(1:10), pve=V1/sum(V1), Method="Standard PCA")
  eval_pcair = mutate(eval_pcair, PC=as.factor(1:10), pve=V1/sum(V1), Method="PC-AIR")
  eval_combined = rbind(eval_pca, eval_pcair)
  
  #Save
  save(eval_combined, file = scree_results)
  
  #Plot
  p_scree <- ggplot(eval_combined, aes(x=PC, y=pve, group=Method, color=Method)) + 
    geom_point() + 
    geom_line() + 
    ylab("Proportion of Variance Explained") + 
    theme_bw(12) +
    ggtitle(paste0("Variance Explained: ", population))
  
  ggsave(file.path(pop_dir, paste0(population, "_scree_plot.png")), p_scree, width = 8, height = 6)
  
  #FINAL SUMMARY
  summary_results <- paste0(results_base, "_summary.RData")
  load(summary_results)
  
  if (!exists("permanova_data")) {
    pca_coords <- all_pca[, c("PC1", "PC2", "PC3", "PC4", "PC5")]
    pcair_coords <- all_pcair[, c("PC1", "PC2", "PC3", "PC4", "PC5")]
    
    permanova_data <- rbind(
      data.frame(
        Method = "Standard PCA",
        SuperPop = all_pca$SuperPop,
        as.data.frame(pca_coords)
      ),
      data.frame(
        Method = "PC-AIR",
        SuperPop = all_pcair$SuperPop,
        as.data.frame(pcair_coords)
      )
    )
  }
  
  #Calculate summary statistics by SuperPop
  summary_stats <- permanova_data %>%
    group_by(Method, SuperPop) %>%
    summarize(
      PC1_mean = mean(PC1),
      PC1_sd = sd(PC1),
      PC2_mean = mean(PC2),
      PC2_sd = sd(PC2),
      PC3_mean = mean(PC3),
      PC3_sd = sd(PC3),
      PC4_mean = mean(PC4),
      PC4_sd = sd(PC4),
      PC5_mean = mean(PC5),
      PC5_sd = sd(PC5)
    )
  
  #Save summary statistics
  save(summary_stats, permanova_data, file = summary_results)
  write.csv(summary_stats, file.path(pop_dir, paste0(population, "_pc_summary_by_superpop.csv")), row.names = FALSE)
  
  #Create density plots
  p_density1 <- ggplot(permanova_data, aes(x=PC1, fill=Method)) + 
    geom_density(alpha=0.5) + 
    theme_bw(12) +
    ggtitle(paste0(population, " - Distribution of PC1 by Method"))
  
  p_density2 <- ggplot(permanova_data, aes(x=PC2, fill=Method)) + 
    geom_density(alpha=0.5) + 
    theme_bw(12) +
    ggtitle(paste0(population, " - Distribution of PC2 by Method"))
  
  #Boxplots
  p_boxplot <- ggplot(permanova_data, aes(x=SuperPop, y=PC1, fill=Method)) + 
    geom_boxplot() + 
    theme_bw(12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste0(population, " - PC1 Distribution by SuperPop and Method"))
  
  #Save
  ggsave(file.path(pop_dir, paste0(population, "_pc1_density.png")), p_density1, width = 8, height = 6)
  ggsave(file.path(pop_dir, paste0(population, "_pc2_density.png")), p_density2, width = 8, height = 6)
  ggsave(file.path(pop_dir, paste0(population, "_pc1_boxplot.png")), p_boxplot, width = 10, height = 6)
  
  return(pop_dir)
}

#Loop through each population and run analysis
for (population in populations) {
  run_population_analysis(population)
}

#Create a directory for cross-population comparisons
cross_pop_dir <- file.path(pcair_dir, "cross_population_comparison")
dir.create(cross_pop_dir, recursive = TRUE)

#Collect PVE data across all populations
all_pve_data <- data.frame()

for (population in populations) {
  scree_results <- file.path(pcair_dir, population, paste0(population, "_analysis_results_scree.RData"))
  load(scree_results)
  eval_combined$Population <- population
  all_pve_data <- rbind(all_pve_data, eval_combined)
}

#Comparison plot
if (nrow(all_pve_data) > 0) {
  # Compare PVE for PC1 across populations
  p_pve_comparison <- ggplot(all_pve_data[all_pve_data$PC == 1,], 
                             aes(x=Population, y=pve, fill=Method)) + 
    geom_bar(stat="identity", position="dodge") + 
    theme_bw(12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Proportion of Variance Explained by PC1") +
    ylab("Proportion of Variance")
  
  ggsave(file.path(cross_pop_dir, "pop_pve_comparison.png"), 
         p_pve_comparison, width = 10, height = 6)
}

#Collect Mantel test results across populations
mantel_results <- data.frame(Population = character(),
                             Correlation = numeric(),
                             Significance = numeric(),
                             stringsAsFactors = FALSE)

for (population in populations) {
  distance_results <- file.path(pcair_dir, population, paste0(population, "_analysis_results_distance.RData"))
  
  tmp_env <- new.env()
  load(distance_results, envir = tmp_env)
  
  if (!is.null(tmp_env$dist_correlation)) {
    mantel_results <- rbind(mantel_results, 
                            data.frame(Population = population,
                                       Correlation = tmp_env$dist_correlation$statistic,
                                       Significance = tmp_env$dist_correlation$signif,
                                       stringsAsFactors = FALSE))
  }
}

if (nrow(mantel_results) > 0) {
  write.csv(mantel_results, file.path(cross_pop_dir, "mantel_comparison.csv"), row.names = FALSE)
  
  #Plot Mantel test correlations across populations
  if (nrow(mantel_results) > 1) {
    p_mantel <- ggplot(mantel_results, aes(x=Population, y=Correlation, fill=Significance < 0.05)) + 
      geom_bar(stat="identity") + 
      theme_bw(12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle("Mantel Test Correlations Between PCA and PC-AIR Distance Matrices") +
      ylab("Correlation") +
      scale_fill_discrete(name="Significant (p<0.05)")
    
    ggsave(file.path(cross_pop_dir, "mantel_correlation_comparison.png"), p_mantel, width = 10, height = 6)
  }
}

#Compare PERMANOVA results across populations
permanova_method_results <- data.frame(Population = character(),
                                       R2 = numeric(),
                                       P_value = numeric(),
                                       stringsAsFactors = FALSE)

for (population in populations) {
  permanova_results <- file.path(pcair_dir, population, paste0(population, "_analysis_results_permanova.RData"))
  
  tmp_env <- new.env()
  load(permanova_results, envir = tmp_env)
    
  if (exists("perm_method", envir = tmp_env)) {
    method_effect <- tmp_env$perm_method
      
    if ("R2" %in% colnames(method_effect) && "Pr(>F)" %in% colnames(method_effect)) {
      permanova_method_results <- rbind(permanova_method_results, 
                                        data.frame(Population = population,
                                                    R2 = method_effect["Method", "R2"],
                                                    P_value = method_effect["Method", "Pr(>F)"],
                                                    stringsAsFactors = FALSE))
    }
  }
}

if (nrow(permanova_method_results) > 0) {
  write.csv(permanova_method_results, file.path(cross_pop_dir, "permanova_method_comparison.csv"), row.names = FALSE)
  
  #Plot PERMANOVA R2 values across populations
  if (nrow(permanova_method_results) > 1) {
    p_permanova <- ggplot(permanova_method_results, aes(x=Population, y=R2, fill=P_value < 0.05)) + 
      geom_bar(stat="identity") + 
      theme_bw(12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle("PERMANOVA R² Values for Method Effect (PCA vs PC-AIR)") +
      ylab("R²") +
      scale_fill_discrete(name="Significant (p<0.05)")
    
    ggsave(file.path(cross_pop_dir, "permanova_r2_comparison.png"), p_permanova, width = 10, height = 6)
  }
}

#Compare PC1 distributions across populations
all_pc1_data <- data.frame()

for (population in populations) {
  summary_results <- file.path(pcair_dir, population, paste0(population, "_analysis_results_summary.RData"))
  
  tmp_env <- new.env()
  load(summary_results, envir = tmp_env)
    
  if (exists("permanova_data", envir = tmp_env)) {
    pc1_data <- tmp_env$permanova_data[, c("Method", "SuperPop", "PC1")]
    pc1_data$Population <- population
    all_pc1_data <- rbind(all_pc1_data, pc1_data)
  }
}

if (nrow(all_pc1_data) > 0) {
  #Create violin plot of PC1 by method and population
  p_violin <- ggplot(all_pc1_data, aes(x=Population, y=PC1, fill=Method)) + 
    geom_violin(position=position_dodge(width=0.9)) + 
    theme_bw(12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("PC1 Distributions by Population and Method")
  
  ggsave(file.path(cross_pop_dir, "pc1_distribution_comparison.png"), 
         p_violin, width = 12, height = 6)
  
  #Create faceted plot by SuperPop
  p_facet <- ggplot(all_pc1_data, aes(x=Method, y=PC1, fill=Method)) + 
    geom_boxplot() + 
    facet_wrap(~SuperPop) +
    theme_bw(12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("PC1 Distributions by SuperPop and Method")
  
  ggsave(file.path(cross_pop_dir, "pc1_by_superpop_faceted.png"), 
         p_facet, width = 14, height = 8)
}