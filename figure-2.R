#
#  Code to reproduce plots in Figure 2
#

source("common.R")


# Read data
hmec_peak_list <- read_hmec_peak_list()

# 13C corrected data for computing 13C enrichment
hmec_mi_data_uncorrected <- readRDS(file.path(mi_data_path, 'hmec_mi_data_censored.rds'))
# uncorrected data
hmec_mi_data <- readRDS(file.path(mi_data_path, 'hmec_mi_data.rds'))

n_peaks <- length(hmec_mi_data$peak_ids)
n_experiments <- length(hmec_mi_data$experiments)

# verify that MIData peak IDs matches peak list
stopifnot(all(hmec_mi_data$peak_ids == hmec_peak_list$peak_id))

# precomputed distance matrix
hmec_dm <- readRDS(file.path(mid_distance_path, 'hmec_13c_corr_dm.rds'))

#
#  Fig 2b 13C enrichment
#

# 13C enrichment values, peaks x experiments

enrichment_mat <- t(
    sapply(1:n_peaks,
        function(p) {
            apply(
               get_avg_mid(hmec_mi_data_uncorrected, p, ), 2,
               isotopic_enrichment
            )
        }
    )
)
dimnames(enrichment_mat) <- list(
    hmec_mi_data$peak_ids,
    hmec_mi_data$experiments
)

# number of experiments with labeling for each peak
table(rowSums(enrichment_mat > 0.05))

# number of metabolite with labeling in > 1 experiment
sum(rowSums(enrichment_mat > 0.05) > 1)


# cluster and plot enrichment matrix
enrichment_clust <- hclust(
    as.dist(1 - cor(t(enrichment_mat))),
    method = "average"
)

image(
    t(enrichment_mat[enrichment_clust$order, ]),
    col = colorRampPalette(c("white", "blue"))(100)
)

#
# ED Fig 2x relative errors across replicates
#

matrix_row_sd <- function(mat) apply(mat, 1, sd)

matrix_row_mean <- function(mat) apply(mat, 1, mean)

matrix_rel_err <- function(mat)
{
   row_means <- matrix_row_mean(mat)
   row_sds <- matrix_row_sd(mat)
   ifelse(row_means > 0, row_sds / row_means, NA)
}

mean_matrix <- sapply(
   hmec_mi_data_uncorrected$experiments,
   function(exp) {
      sapply(
         1:n_peaks,
         function(peak) {
            matrix_row_mean(get_mids(hmec_mi_data_uncorrected, peak, exp))
         }
      ) %>% unlist()
   }
)

rel_err_matrix <- sapply(
   hmec_mi_data_uncorrected$experiments,
   function(exp) {
      sapply(
         1:n_peaks,
         function(peak) {
            matrix_rel_err(get_mids(hmec_mi_data_uncorrected, peak, exp))
         }
      ) %>% unlist()
   }
)

# relative errors for "detectable" MIs (mean fraction > 5%)
hist(pmin(rel_err_matrix[mean_matrix > 0.05], 0.3), n = 100)

# fraction of MIs with relative error < 5%
sum(rel_err_matrix[mean_matrix > 0.05] < 0.05, na.rm = TRUE) / sum(mean_matrix > 0.05, na.rm = TRUE)


#
# Fig 2c glutathione
#

# selected experiments to plot
selected_exp <- c("cys", "glc", "gln", "glu", "gly", "ser")

# 5108 glutathione +H
mids_gthrd <- get_mid_matrix(hmec_mi_data, "5108", selected_exp)
plot_mid_matrix(mids_gthrd, max_mi_fraction = 0.3)


#
# Fig 2d inosine
#

# selected experiments to plot
selected_exp <- c("glc", "gly", "ser")

# 6795 hypoxantine
mids_hxan <- get_mid_matrix(hmec_mi_data, "6795", selected_exp)
plot_mid_matrix(mids_hxan, max_mi_fraction = 0.3)

# 1679 ribose-5P (pentose)
mids_pentose <- get_mid_matrix(hmec_mi_data, "1679", selected_exp)
plot_mid_matrix(mids_pentose, max_mi_fraction = 0.3)

# 1467 inosine
mids_ins <- get_mid_matrix(hmec_mi_data, "1467", selected_exp)
plot_mid_matrix(mids_ins, max_mi_fraction = 0.3)

# convolution hxan x pentose
mids_conv <- convolute_all(mids_hxan, mids_pentose)
plot_mid_matrix(mids_conv, max_mi_fraction = 0.3)


#
# ED Figure 2b MID distance distribution
#

# histogram of distances, NA set to maximum
hist(hmec_dm[lower.tri(hmec_dm)], n = 100)


#
# Figure 2e UMap projection
#

umap_proj <- umap_projection(
   hmec_dm, n_neighbors = 15, random_seed = 571632932)

fig2e_clusters <- read_tsv(file.path(input_data_path, "fig_2_clusters.tsv")) %>%
    mutate(peak_id = as.character(peak_id)) %>%
    mutate(
        color = case_when(
            cluster_name %in% c("asn", "glu/gln", "ser lipids", "sugars", "val") ~"#0003FF",
            cluster_name %in% c("glutathione", "sulfur amino acids", "trp (indole)") ~"#50F200",
            cluster_name %in% c("leu/ile", "purines", "pyrimidines", "thr") ~"#7854A2",
            cluster_name %in% c("his", "lysine", "TCA & asp") ~"#DF1B53",
            cluster_name %in% c("gln", "ser") ~"#ABBDFF"))

fig2e_colors <- hmec_peak_list %>%
    left_join(fig2e_clusters, by = 'peak_id') %>%
    mutate(color = coalesce(color, "#BFBFBF")) %>%
    pull(color)

plot_umap(umap_proj, colors = fig2e_colors)



# export peak list with highlighted clusters and UMap coordinates, for Suppl Table 1
hmec_peak_list %>%
   select(peak_id) %>%
   inner_join(umap_proj, by = 'peak_id') %>%
   left_join(
      fig2e_clusters %>% select(peak_id, cluster_name),
      by = 'peak_id') %>%
   write.table(
      file.path('suppl_table_1_fig2e_columns.tsv'),
      sep = '\t', row.names = FALSE, quote = FALSE
   )

# interactive plot with peak numbers

plotly_tooltips <- read_plotly_tooltips()
plot_umap_interactive(umap_proj, plotly_tooltips$tooltip)


#
# ED Figure 2c  UMap distance vs MID distance
#

# all pairwise distances, UMap vs MID distance. this generally does not agree well
# pairwise euclidean distances of 2D UMap projection
umap_dist <- umap_proj %>% select(umap_1, umap_2) %>% dist()
local({
    index <- sample(length(umap_dist), 10000)
    plot(as.dist(hmec_dm)[index], umap_dist[index], col = alpha("black", 0.4))
})

# Spearman correlation (lower diagonal)
cor(as.dist(hmec_dm), umap_dist, method = "spearman")


#
# Figure 2g citryl-glutamate
#

# selected experiments to plot
selected_exp <- c("glc", "gln", "glu")

# citrate
mids_cit <- get_mid_matrix(hmec_mi_data, "1917", selected_exp)
plot_mid_matrix(c13correct_cols(mids_cit), max_mi_fraction = 0.3)

# glutamate
mids_glu <- get_mid_matrix(hmec_mi_data, "6683", selected_exp)
plot_mid_matrix(c13correct_cols(mids_glu), max_mi_fraction = 0.3)

# 1171 citryl-glutamate (candidate)
mids_citglu <- get_mid_matrix(hmec_mi_data, "1171", selected_exp)
plot_mid_matrix(c13correct_cols(mids_citglu), max_mi_fraction = 0.3)

# convolution of citrate and glutamate
plot_mid_matrix(
    c13correct_cols(convolute_all(mids_cit, mids_glu)),
    max_mi_fraction = 0.3)

#
# Figure 2h MID-based network
#

dist_cutoff <- 0.3
pair_index <- which(hmec_dm < 0.3, arr.ind = TRUE)

# table of node pairs
network_table <- as.data.frame(pair_index) %>%
   rename(index_1 = row, index_2 = col) %>%
   mutate(
      peak_id_1 = hmec_mi_data$peak_ids[index_1],
      peak_id_2 = hmec_mi_data$peak_ids[index_2],
      distance = hmec_dm[pair_index],
      n_carbon_1 = hmec_mi_data$peak_n_atoms[index_1],
      n_carbon_2 = hmec_mi_data$peak_n_atoms[index_2]
   ) %>%
   filter(index_1 < index_2)

# number of metabolites with at least one connection
network_nodes <- union(network_table$peak_id_1, network_table$peak_id_2)
network_nodes %>% length

# average node degree in the network
nrow(network_table) / length(network_nodes)
    
# export table for visualization in Cytoscape
write.table(
   network_table,
   "all-tracers-13c-corr-dist0.3-edges.tsv",
   col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE
)

