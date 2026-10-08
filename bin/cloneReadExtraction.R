#!/usr/bin/env Rscript

# This script accepts a .txt file of somatic hypermutation SNVs with VAF annotations and
# a .tsv of per-base information of all reads containing SHM SNVs. It outputs a diagnostic
# plot of BMix binomial fit to the distribution of VAFs and .txt files for each set of unique
# reads containing SHM SNVs at each VAF cluster.
# author: pblaney

#########################
#####   Libraries   #####

suppressPackageStartupMessages(library(BMix))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

set.seed(999)

#########################
#####   Execution   #####

# Accept command line arguments as input
input_args <- commandArgs(trailingOnly = T)

somatic_hypermutation_snvs_txt <- input_args[1]

reads_tsv <- input_args[2]

sample_id <- stringr::str_remove(string = basename(somatic_hypermutation_snvs_txt),
                                 pattern = "\\..*\\.txt")

cat("\n")
cli::cli_alert_info("Beginning {crayon::cyan('BMix')} workflow run for {crayon::red({sample_id})}")

# Read in SHM SNV data and filter out any mutations with a VAF less than 0.05 to
# remove hyper sub-clonal mutations non-likely to capture 
shm_muts <- data.table::fread(file = somatic_hypermutation_snvs_txt,
                              sep = "\t",
                              header = FALSE,
                              col.names = c("chrom","pos","ref","alt","sample","AD","DP","VAF")) %>%
                          dplyr::filter(VAF > 0.05)

# Build data.frame for BMix clustering
# BMix expects whole integer values so need to convert VAF to value between 0-100 instead of 0-1
bmix_vafs <- data.table::data.table("vaf" = as.integer(shm_muts$VAF * 100),
                                    "trials" = 100)

# Run BMix to identify clusters in density of VAFs
# Max of 4 possible VAF clusters in binomial fit
bmix_vaf_clusters <- BMix::bmixfit(data = bmix_vafs,
                                   K.Binomials = 1:4,
                                   K.BetaBinomials = 0)

# Isolate the cluster VAF means and rank them in ascending order
ranked_vaf_cluster_means <- sort(bmix_vaf_clusters$B.params)

# Set colors for plotting based on number of clusters
cluster_color_pal <- c("red","blue","green","purple")

# Plot the clusters along the density curve
cat("\n")
cli::cli_alert_info("Generating diagnostic plot")
cat("\n")
bmix_fit_plot <- BMix::plot_density(x = bmix_vaf_clusters, data = bmix_vafs) +
  ggplot2::scale_color_manual(values = cluster_color_pal[1:length(ranked_vaf_cluster_means)]) +
  ggplot2::theme(panel.background = element_blank(),
                 panel.border = element_rect(fill = "transparent"),
                 panel.grid.major.x = element_blank(),
                 panel.grid.major.y = element_line(color = "gray60", linetype = "dashed"),
                 panel.grid.minor.x = element_blank(),
                 panel.grid.minor.y = element_blank(),
                 legend.position = "inside",
                 legend.position.inside = c(0.88,0.85),
                 legend.text = element_text(size = 9),
                 legend.background = element_rect(fill = alpha("gray80", alpha = 0.8)))

# combine into single DT the mutations, position , VAF cluster
clutered_shm_muts <- shm_muts %>%
                       dplyr::select(chrom,pos,alt,VAF) %>%
                       dplyr::mutate("cluster" = bmix_vaf_clusters$labels,
                                     "mut_key" = paste(chrom, pos, alt, sep = ":"))

# Plot the location of each mutation with VAF and cluster ID
clustered_mut_by_pos <- ggplot2::ggplot(clutered_shm_muts) +
  ggplot2::geom_point(aes(x = pos, y = VAF, colour = cluster), size = 1.5) +
  ggplot2::scale_color_manual(name = "Cluster", values = cluster_color_pal[1:length(ranked_vaf_cluster_means)]) +
  ggplot2::labs(x = paste0("Position on ", stringr::str_to_title(unique(clutered_shm_muts$chrom)))) +
  ggplot2::theme(panel.background = element_blank(),
                 panel.border = element_rect(fill = "transparent"),
                 panel.grid.major.x = element_blank(),
                 panel.grid.major.y = element_line(color = "gray40", linetype = "dotted"),
                 panel.grid.minor.x = element_blank(),
                 panel.grid.minor.y = element_blank(),
                 legend.position = "inside",
                 legend.position.inside = c(0.12,0.80),
                 legend.text = element_text(size = 8),
                 legend.background = element_rect(fill = alpha("gray80", alpha = 0.8)))

# diagnostic plot
bmix_diagnostic_plot <- bmix_fit_plot + clustered_mut_by_pos

# Save the BMix plot
ggplot2::ggsave(filename = paste0(sample_id, "_bmix_fit_diagnostic.pdf"),
                plot = bmix_diagnostic_plot,
                path = getwd(),
                width = 250,
                height = 125,
                units = "mm",
                device = "pdf")

# Read in sam2tsv file with read names
cat("\n")
cli::cli_alert_info("Extracting reads per VAF cluster")
shm_reads_by_base <- data.table::fread(file = reads_tsv,
                                       sep = "\t",
                                       header = TRUE) %>%
                       dplyr::mutate("mut_key" = paste(CHROM,`REF-POS1`,`READ-BASE`,sep = ":"))

# Gather all reads with mutations per cluster
reads_per_clustered_mut <- data.table::data.table()
for(i in 1:dplyr::n_distinct(clutered_shm_muts$cluster)) {
  
  # Isolate all mutations per cluster
  per_cluster_shm_muts <- clutered_shm_muts %>%
                            dplyr::filter(cluster == unique(clutered_shm_muts$cluster)[i])
  
  # Find all reads that have the mutation present
  for(j in 1:nrow(per_cluster_shm_muts)) {
    reads_per_clustered_mut <- rbind(reads_per_clustered_mut,
                                        shm_reads_by_base %>%
                                            dplyr::filter(mut_key == per_cluster_shm_muts$mut_key[j]) %>%
                                            dplyr::mutate("cluster" = unique(per_cluster_shm_muts$cluster)) %>%
                                            dplyr::select(`#Read-Name`,mut_key,cluster))
  }
}

# Need to find reads exclusive to each cluster
# To achieve this will start with most clonal, then second, third, etc.
# For each subsequent cluster, the reads must not be defined in the others
list_of_cluster_reads <- list()
for(i in 1:dplyr::n_distinct(reads_per_clustered_mut$cluster)) {
   cluster_unique_reads <- reads_per_clustered_mut %>%
                             dplyr::filter(cluster == names(ranked_vaf_cluster_means[i])) %>%
                             dplyr::distinct(`#Read-Name`)

   # If the first iteration, which corresponds to the most sub-clonal cluster, keep all reads
   # For all other subsequent clusters, the reads must not be present in previous clusters
   if(i == 1) {
     list_of_cluster_reads[[i]] <- cluster_unique_reads
   } else {
     list_of_cluster_reads[[i]] <- cluster_unique_reads %>%
                                     dplyr::filter(!`#Read-Name` %in% data.table::rbindlist(list_of_cluster_reads[-i])$`#Read-Name`)
  }
}

# Check how many reads are left in each cluster,
# only export them if greater than 5
cli::cli_alert_info("Writing file of read QNAMEs")
for(i in 1:length(list_of_cluster_reads)) {
  if(nrow(list_of_cluster_reads[[i]]) > 5) {
    data.table::fwrite(x = list_of_cluster_reads[[i]],
                       file = paste0(sample_id, "_",
                                     stringr::str_replace(string = names(ranked_vaf_cluster_means[i]), pattern = " ", replacement = ""), "_",
                                     round(ranked_vaf_cluster_means[i], digits = 3), "_read_qnames.txt"),
                       sep = "\t",
                       col.names = FALSE)
  } else if(nrow(list_of_cluster_reads[[i]]) <= 5) {
    cli::cli_alert_info("NOTE: {names(ranked_vaf_cluster_means[i])} was excluded as it had less than 5 unique reads")
  }
}

cli::cli_alert_success("{crayon::green('Step complete ...')}")
cat("\n")
