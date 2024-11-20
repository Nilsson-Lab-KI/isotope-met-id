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

# subset MI data to TCA metabolites and compute specific distance matrix
tca_met_ids <- c("cit_m", "icit_m", "akg_m", "succoa_m", "succ_m", "fum_m", "mal-L_m", "oaa_m", "accoa_m", "co2_m")
sim_mi_data <- readRDS(sim_mi_data_path(mi_stdev, 1))
sim_mi_data_tca <- midata_subset(sim_mi_data, tca_met_ids)

stopifnot(all(sim_mi_data_tca$peak_ids == tca_met_ids))

assign_list[sim_dm_tca, sim_conv_index_tca] <- conv_reduce_all(
   sim_mi_data_tca,
   1:length(sim_mi_data_tca$experiments),
   f = midist::euclidean_sum_dist,
   g = which.min
)
sim_dm_tca <- replace_na_with_max(sim_dm_tca)

# exclude accoa and co2 from the UMap (small convolutants)
tca_met_ids_to_plot <- c("cit_m", "icit_m", "akg_m", "succoa_m", "succ_m", "fum_m", "mal-L_m", "oaa_m")

umap_proj <- umap_projection(
   sim_dm_tca[tca_met_ids_to_plot, tca_met_ids_to_plot],
   n_neighbors = 5, random_seed = 381691)

plot_umap(umap_proj)

plot_umap_interactive(umap_proj, tca_met_ids_to_plot)


#
# Fig 3e precision-recall curves for various noise levels
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
# Figure 3f precision-recall curves for subsets of metabolites
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
# Figure 3g Area under PR curve (AUPR) for pairs, individual tracers
#

gold_standard <- readRDS(file.path(gold_standard_path, 'gold_standard.rds'))

mi_stdev <- 0.01
n_replicates <- 10
sim_dm <- readRDS(sim_dm_ind_exp_path(mi_stdev, "glc", 1))


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

gold_standard <- readRDS(file.path(gold_standard_path, 'gold_standard.rds'))

ranks_experiment <- function(experiment)
{
   sim_dm <- readRDS(sim_dm_ind_exp_path(stdev = 0.01, experiment, rep_nr = 1))
   sim_dm_ranks <- distance_ranks(sim_dm)
   median_ranks(sim_dm_ranks, gold_standard)
}

# this yields a metabolites x experiments matrix
# contains NAs for a few metabolites that have no neighbors in gold standard
rank_threshold <- 20
ranks_matrix <- sapply(c(experiments, "all"), ranks_experiment) %>% replace_na_with_max()
ranks_matrix <- pmin(ranks_matrix, rank_threshold+1)

# remove "undiscoverable" metabolites
ranks_matrix <- ranks_matrix[
   apply(ranks_matrix, MARGIN = 1, min) <= rank_threshold, ]

# cluster metabolites
metabolite_clust <- ranks_matrix[, 1:n_experiments] %>% dist() %>% hclust(method = "ward.D")
metabolite_order <- metabolite_clust$order

# cluster individual tracers, keep "all" in last column
experiment_clust <- ranks_matrix[, 1:n_experiments] %>% t() %>%
   dist() %>% hclust(method = "ward.D")
experiment_order <-  c(experiment_clust$order, n_experiments+1)

plot(experiment_clust)

ranks_matrix[metabolite_order, rev(experiment_order)] %>%
   image(col = colorRampPalette(c("red", "white"))(rank_threshold + 1))


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
# ED Figure 3d Successively adding tracing experiments
#

# here we need to compute AUPR of n matrices in step 1,
# n-1 matrices in step 2, ... 1 matrix in step n, total n*(n+1)/2

sim_mi_data <- readRDS(sim_mi_data_path(stdev = 0.01, rep_nr = 1))

next_best_experiment <- function(experiments, prev_experiments)
{
   try_experiments <- setdiff(experiments, prev_experiments)
   print(try_experiments)
   best_aupr <- 0
   for(experiment in try_experiments) {
      new_experiments <- c(prev_experiments, experiment)
      cat(paste(new_experiments, sep = "+"))
      assign_list[dm, conv_index] <- conv_reduce_all(
         sim_mi_data,
         e = new_experiments,
         f = ifelse(length(new_experiments) > 1, euclidean_sum_dist, euclidean_dist),
         g = which.min
      )
      new_aupr <- continuous_accuracy(dm, gold_standard) %>% aupr()
      cat(paste0(" AUPR = ", round(new_aupr, 3), "\n"))
      if(new_aupr > best_aupr) {
         best_experiment <- experiment
         best_aupr <- new_aupr
      }
   }
   return(data.frame(experiment = best_experiment, aupr = best_aupr))
}

{
   best_experiments <- data.frame(
      iteration = integer(),
      experiment = character(),
      aupr = numeric()
   )
   
   for(i in 1:n_experiments) {
      best_experiments <- bind_rows(
         best_experiments,
         next_best_experiment(
            experiments = experiments,
            prev_experiments = best_experiments$experiment) %>% mutate(iteration = i)
      )
   }
}

best_experiments %>%
   mutate(experiment = factor(experiment, levels = experiment)) %>%
   ggplot(aes(x = experiment, y = aupr)) + 
   geom_point() +
   ylim(0, 0.5) + theme_classic()

