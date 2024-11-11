#
#  Figure 3

source("common.R")

# library(reshape2)


#
# Figure 3c
#

mi_stdev <- 0.01
sim_dm <- readRDS(sim_dm_path(mi_stdev, 1))

umap_proj <- umap_projection(
   sim_dm, n_neighbors = 15, random_seed = 3731594)

# TODO: provide table of cluster annotations 

# fig3c_clusters <- read_tsv(file.path(input_data_path, "fig_2_clusters.tsv")) %>%
#    mutate(peak_id = as.character(peak_id)) %>%
#    mutate(
#       color = case_when(
#          cluster_name %in% c("asn", "glu/gln", "ser lipids", "sugars", "val") ~"#0003FF",
#          cluster_name %in% c("glutathione", "sulfur amino acids", "trp (indole)") ~"#50F200",
#          cluster_name %in% c("leu/ile", "purines", "pyrimidines", "thr") ~"#7854A2",
#          cluster_name %in% c("his", "lysine", "TCA & asp") ~"#DF1B53",
#          cluster_name %in% c("gln", "ser") ~"#ABBDFF"))

# fig2e_colors <- hmec_peak_list %>%
#    left_join(fig2e_clusters, by = 'peak_id') %>%
#    mutate(color = coalesce(color, "#BFBFBF")) %>%
#    pull(color)

plot_umap(umap_proj)

# interactive plot with peak numbers
#plotly_tooltips <- read_plotly_tooltips()
plot_umap_interactive(umap_proj, rownames(sim_dm))

#
# Figure 3d TCA cycle
#

# version w/o accoa and co2 looks quite good

tca_met_ids <- c("cit_m", "icit_m", "akg_m", "succoa_m", "succ_m", "fum_m", "mal-L_m", "oaa_m", "accoa_m", "co2_m")

# subset MI data to TCA metabolites and compute specific distance matrix
sim_mi_data <- readRDS(sim_mi_data_path(mi_stdev, 1))
sim_mi_data_tca <- midata_subset(sim_mi_data, tca_met_ids)

assign_list[sim_dm_tca, sim_conv_index_tca] <- conv_reduce_all(
   sim_mi_data_tca,
   1:length(sim_mi_data_tca$experiments),
   f = midist::euclidean_sum_dist,
   g = which.min
)

# impute missing values with maximal distance
sim_dm_tca[which(is.na(sim_dm_tca))] <- max(sim_dm_tca, na.rm = TRUE)

umap_proj <- umap_projection(
   sim_dm_tca,
   n_neighbors = 4, random_seed = 3629136)

plot_umap(umap_proj)

plot_umap_interactive(umap_proj, tca_met_ids)

# explore MIDS

plot_tca_mid <- function(met_id)
{
   plot_mid_matrix(
      get_mid_matrix(sim_mi_data, met_id, sim_mi_data_tca$experiments),
      max_mi_fraction = 0.3, plot_title = met_id
   )
}

plot_tca_mid("cit_m")
plot_tca_mid("icit_m")
plot_tca_mid("oaa_m")

#
# Figure 3e fatty acid synthesis
#


fas_met_ids <- c("acACP_c", "btACP_c", "hexACP_c", "ocACP_c", "dcaACP_c", "ddcaACP_c", "myrsACP_c", "palmACP_c") # "malcoa_c"
fas_tooltips <- c("C2", "C4", "C6", "C8", "C10", "C12", "C14", "C16")

umap_proj <- umap_projection(
   sim_dm[fas_met_ids, fas_met_ids],
   n_neighbors = 6, random_seed = 3629136)

plot_umap(umap_proj)

plot_umap_interactive(umap_proj, fas_tooltips)


#
# Fig 3f precision-recall curves for various noise levels
#

gold_standard <- readRDS(file.path(gold_standard_path, 'gold_standard.rds'))

mi_stdev_range <- c(0.01, 0.05, 0.1)
n_replicates <- 10
downsample_factor <- 20    # reduce number of pairs plotted by this factor

all(colnames(sim_dm) == colnames(gold_standard))

accuracy_rep <- function(mi_stdev, rep_nr)
{
   sim_dm <- readRDS(sim_dm_path(mi_stdev, rep_nr))
   continuous_accuracy(sim_dm, gold_standard) %>%
      tibble::rowid_to_column('pair_rank') %>%
      filter((pair_rank - 1) %% downsample_factor == 0) %>%
      mutate(mi_stdev = mi_stdev) %>%
      mutate(rep_nr = rep_nr)
}

accuracy_stdev <- function(mi_stdev)
{
   bind_rows(
      lapply(
         1:n_replicates,
         function(rep_nr) accuracy_rep(mi_stdev, rep_nr)))
}

accuracy <- bind_rows(
   lapply(mi_stdev_range, accuracy_stdev),
   # add noise-free accuracy
   accuracy_rep(0, 1)
)

accuracy %>% filter(mi_stdev == 0) %>% filter(rep_nr == 1) %>% head()

# average precision and recall values over replicate samples
mean_accuracy <- accuracy %>%
   group_by(mi_stdev, pair_rank) %>%
   summarize(
      mean_precision = mean(precision),
      mean_recall = mean(recall))

mean_accuracy %>%
   filter(pair_rank > 1) %>%
   ggplot(aes(x = mean_recall, y = mean_precision, group = mi_stdev, color = factor(mi_stdev))) +
      geom_line() +
      theme_classic() + xlim(0, 1) + ylim(0, 1)


#
# Figure 3g precision-recall curves for subsets of metabolites
#

mi_stdev <- 0.01
subset_sizes <- c(300, 200, 100, 50)
n_replicates <- 10
max_curve_points <- 1000

accuracy_subset_rep <- function(subset_size, rep_nr)
{
   subset_dm <- readRDS(subset_dm_path(subset_size, rep_nr))
   subset_gold_standard <- gold_standard[rownames(subset_dm), colnames(subset_dm)]
   
   n_pairs <- nrow(subset_dm) * (nrow(subset_dm) - 1) / 2
   downsample_factor <- ceiling(n_pairs / max_curve_points)
   continuous_accuracy(subset_dm, subset_gold_standard) %>%
      tibble::rowid_to_column('pair_rank') %>%
      filter((pair_rank - 1) %% downsample_factor == 0) %>%
      mutate(subset_size = subset_size) %>%
      mutate(rep_nr = rep_nr)
}

accuracy_subset <- function(subset_size)
{
   bind_rows(
      lapply(
         1:n_replicates,
         function(rep_nr) accuracy_subset_rep(subset_size, rep_nr)))
}

accuracy <- bind_rows(lapply(subset_sizes, accuracy_subset))

# average precision and recall values over replicate samples
mean_accuracy <- accuracy %>%
   group_by(subset_size, pair_rank) %>%
   summarize(
      mean_precision = mean(precision),
      mean_recall = mean(recall))

mean_accuracy %>%
   filter(pair_rank > 1) %>%
   ggplot(aes(x = mean_recall, y = mean_precision, group = subset_size, color = factor(subset_size))) +
   geom_line() +
   theme_classic() + xlim(0, 1) + ylim(0, 1)

#
# ED Figure 3b
#

fraction_derived <- readRDS(file.path(gold_standard_path, 'fraction_derived.rds'))
metabolite_pathway <- read_metabolite_pathway()

frac_derived_long <- melt(
    fraction_derived,
    varnames = c("metabolite_1", "metabolite_2"),
    value.name = "fraction_derived"
)

# pairs of metabolites in the same pathway, lower-triangular
pathway_pairs <- metabolite_pathway %>%
    rename(metabolite_1 = id) %>%
    inner_join(
        metabolite_pathway %>%
            rename(metabolite_2 = id),
        by = 'pathway',
        relationship = 'many-to-many'
    ) %>%
    filter(metabolite_1 < metabolite_2)


# fraction derived between metabolite in same pathway, ignoring compartments
frac_derived_same_pathway <- frac_derived_long %>%
    mutate(
        metabolite_1 = str_sub(metabolite_1, end = -3),
        metabolite_2 = str_sub(metabolite_2, end = -3)
    ) %>%
    inner_join(
        pathway_pairs,
        by = c("metabolite_1", "metabolite_2"),
        relationship = 'many-to-many'
    )

frac_derived_diff_pathway <- frac_derived_long %>%
    anti_join(
        frac_derived_same_pathway,
        by = c("metabolite_1", "metabolite_2")
    )
    
# plot histograms
hist(
    frac_derived_same_pathway$fraction_derived,
    breaks = 50, freq = FALSE)

hist(
    frac_derived_diff_pathway$fraction_derived,
    breaks = 50, freq = FALSE)


#
# ED Figure 3c Comparison with path lengths
#
# for computation of path lengths from atom map, see Mathematica code
#

path_lengths = readRDS(file.path(gold_standard_path, 'path_lengths.rds'))

local({
    original_par = par(no.readonly = TRUE)
    par(pty = "s")
    plot(
        x = fraction_derived,
        y = path_lengths,
        col = alpha("black", 0.05)
    )
    par(original_par)
})

