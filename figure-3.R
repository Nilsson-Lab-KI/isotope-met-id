#
#  Figure 3

source("common.R")

# library(reshape2)

experiments <- readRDS(file.path(mi_data_path, "simulated_mi_data_t75.rds"))$experiments
n_experiments <- length(experiments)

#
# Figure 3c
#

mi_stdev <- 0.01
sim_dm <- readRDS(sim_dm_path(mi_stdev, 1))

umap_proj <- umap_projection(
   sim_dm, n_neighbors = 15, random_seed = 3731594)

# generate color for pathway annotations

n_pathways = metabolite_pathway %>% distinct(pathway) %>% nrow()

metabolite_pathway <- read_metabolite_pathway() %>%
   inner_join(
      data.frame(
         peak_id = umap_proj$peak_id,
         id = str_sub(umap_proj$peak_id, end = -3)),
      by = "id") %>%
   group_by(pathway) %>%
   mutate(pathway_index = cur_group_id()) %>%
   ungroup() %>%
   mutate(color = hsv(pathway_index / n_pathways, 1, 1))

fig_3c_colors <- umap_proj %>%
   select(peak_id) %>%
   left_join(
      metabolite_pathway %>% select(peak_id, pathway, color),
      by = 'peak_id') %>%
   mutate(color = coalesce(color, "#BFBFBF"))

umap_proj %>%
   ggplot(aes(x = umap_1, y = umap_2, colour = fig_3c_colors$pathway)) +
   geom_point(alpha = 0.7, show.legend = TRUE) +
   theme_classic()

plot_umap(umap_proj, colors = fig_3c_colors) + scale_color

# interactive plot with peak numbers
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

# plot_tca_mid <- function(met_id)
# {
#    plot_mid_matrix(
#       get_mid_matrix(sim_mi_data, met_id, sim_mi_data_tca$experiments),
#       max_mi_fraction = 0.3, plot_title = met_id
#    )
# }
# 
# plot_tca_mid("cit_m")
# plot_tca_mid("icit_m")
# plot_tca_mid("oaa_m")

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

accuracy_rep <- function(mi_stdev, rep_nr)
{
   sim_dm <- readRDS(sim_dm_path(mi_stdev, rep_nr))
   continuous_accuracy(sim_dm, gold_standard) %>%
      tibble::rowid_to_column('pair_rank') %>%
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

# average precision and recall values over replicate samples
mean_accuracy <- accuracy %>%
   group_by(mi_stdev, pair_rank) %>%
   summarize(
      mean_precision = mean(precision),
      mean_recall = mean(recall))

# plot curves, simplified to ~100 points
mean_accuracy %>%
   mutate(mean_recall = round(mean_recall, digits = 2)) %>%
   distinct(mean_recall, .keep_all = TRUE) %>%
   ggplot(aes(x = mean_recall, y = mean_precision, group = mi_stdev, color = factor(mi_stdev))) +
      geom_line() +
      theme_classic() + xlim(0, 1) + ylim(0, 1)

   
#
# Figure 3g precision-recall curves for subsets of metabolites
#

mi_stdev <- 0.01
subset_sizes <- c(300, 200, 100, 50)
n_replicates <- 10

accuracy_subset_rep <- function(subset_size, rep_nr)
{
   subset_dm <- readRDS(subset_dm_path(subset_size, rep_nr))
   subset_gold_standard <- gold_standard[rownames(subset_dm), colnames(subset_dm)]
   
   continuous_accuracy(subset_dm, subset_gold_standard) %>%
      tibble::rowid_to_column('pair_rank') %>%
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
   mutate(mean_recall = round(mean_recall, digits = 2)) %>%
   distinct(mean_recall, .keep_all = TRUE) %>%
   ggplot(aes(x = mean_recall, y = mean_precision, group = subset_size, color = factor(subset_size))) +
   geom_line() +
   theme_classic() + xlim(0, 1) + ylim(0, 1)


#
# Figure 3h Area under PR curve (AUPR) for pairs, individual tracers
#

gold_standard <- readRDS(file.path(gold_standard_path, 'gold_standard.rds'))

mi_stdev <- 0.01
n_replicates <- 10
sim_dm <- readRDS(sim_dm_ind_exp_path(mi_stdev, "glc", 1))

experiments <- readRDS(sim_mi_data_path(mi_stdev, 1))$experiments
n_experiments <- length(experiments)

# compute AUPR from an accuracy data.frame with columns 'recall', 'precision'
# rows are assumed to be sorted by the 'recall' column
aupr <- function(accuracy)
{
   n <- nrow(accuracy)
   sum(
      diff(accuracy$recall) * (accuracy$precision[-n] + accuracy$precision[-1]) * 0.5
   )
}

aupr_rep <- function(experiment, rep_nr)
{
   sim_dm <- readRDS(sim_dm_ind_exp_path(mi_stdev, experiment, rep_nr))
   data.frame(
      experiment = experiment,
      rep_nr = rep_nr,
      aupr = continuous_accuracy(sim_dm, gold_standard) %>% aupr()
   )
}

aupr_experiment <- function(experiment)
{
   bind_rows(
      lapply(
         1:n_replicates,
         function(rep_nr) aupr_rep(experiment, rep_nr)))
}

aupr_all <- bind_rows(
   lapply(c(experiments, "all"), aupr_experiment),
)

aupr_all_mean_sd <- aupr_all %>%
   mutate(experiment = factor(experiment, levels = experiments)) %>%
   group_by(experiment) %>%
   summarize(
      aupr_mean = mean(aupr),
      aupr_sd = sd(aupr))

aupr_all_mean_sd %>%
   ggplot(aes(experiment, aupr_mean, fill = experiment)) +
      geom_bar(
         stat = "identity",
         position = position_dodge(width = 0.8), width = 0.7, fill = "gray") +
      geom_errorbar(
         aes(ymin = aupr_mean - aupr_sd, ymax = aupr_mean + aupr_sd),
         position = position_dodge(width = 0.8), width = 0.25) +
      labs(x = "Experiment", y = "AUPR") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))


#
# Figure 3i metabolite ranks (all tracers)
#

mi_stdev <- 0.01
sim_dm <- readRDS(sim_dm_path(mi_stdev, rep_nr = 1))

# matrix where each column is a list of neighbor ranks
distance_ranks <- function(dm)
{
   diag(dm) <- NA
   apply(
      dm, MARGIN = 2,
      function(col) rank(col, na.last = TRUE, ties.method = "average")
   )
}

# median ranks of gold standard neighbors
median_ranks <- function(rank_matrix, gold_standard)
{
   apply(
      rank_matrix * gold_standard,
      MARGIN = 2,
      function(col) median(col[col != 0])
   )
}

permute_columns <- function(mat)
{
   apply(mat, MARGIN = 2, function(col) sample(col, size = nrow(mat)))
}

{
   set.seed(84771582)
   sim_dm_ranks <- distance_ranks(sim_dm)
   metabolite_ranks = data.frame(
      median_rank = median_ranks(sim_dm_ranks, gold_standard),
      random_rank = median_ranks(permute_columns(sim_dm_ranks), gold_standard)
   )
}

metabolite_ranks_long <- metabolite_ranks %>%
   arrange(median_rank) %>%
   tibble::rowid_to_column(var = "index") %>%
   pivot_longer(cols = c('median_rank', 'random_rank'), values_to = 'rank', names_to = 'rank_type')

metabolite_ranks_long %>%
   ggplot(aes(x = index, y = rank, group = rank_type, color = rank_type)) +
      geom_point() +
      labs(x = NULL, y = "Median rank") + 
      scale_y_reverse() +
      scale_color_manual(values = c("red", "black")) +
      theme_classic() 


#
# Figure 3j metabolite ranks per experiment
#

# TEMPORARY: test with Deniz' data

load(
   file.path("C:/tmp/isotope-met-id-deniz", "data", "reverse-engineering",
             "Results", "simulation", "Rdata_files", "t_75_low_noise_sample_1_N404.RData"))
names(matrices) <- c("all", experiments)


ranks_experiment <- function(experiment)
{
   sim_dm <- readRDS(sim_dm_ind_exp_path(stdev = 0.01, experiment, rep_nr = 1))
   #sim_dm <- matrices[[experiment]]$distance_matrix
   sim_dm_ranks <- distance_ranks(sim_dm)
   median_ranks(sim_dm_ranks, gold_standard)
}

# this yields a metabolites x experiments matrix
# contains NAs for a few metabolites that have no neighbors in gold standard
rank_threshold <- 20
ranks_matrix <- sapply(c(experiments, "all"), ranks_experiment) %>% replace_na_with_max()
ranks_matrix <- pmin(ranks_matrix, rank_threshold+1)
#ranks_matrix <- 1 - pmin(ranks_matrix, rank_threshold)/rank_threshold

# cluster metabolites
metabolite_clust <- ranks_matrix[, 1:n_experiments] %>% dist() %>% hclust(method = "ward.D")
metabolite_order <- metabolite_clust$order

# cluster individual tracers, keep "all" in last column
experiment_clust <- ranks_matrix[, 1:n_experiments] %>% t() %>%
   dist() %>% hclust(method = "ward.D")
experiment_order <-  c(experiment_clust$order, n_experiments+1)

plot(experiment_clust)

# TODO: export plot_matrix from midist
ranks_matrix[metabolite_order, rev(experiment_order)] %>% t() %>%
   midist:::plot_matrix() +
      theme(axis.ticks = element_blank(), axis.text.x = element_blank())

image(
   ranks_matrix[rev(metabolite_order), rev(experiment_order)],
   col = colorRampPalette(c("red", "white"))(rank_threshold)
)

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

