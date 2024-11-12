#
# Compute MID distances matrices from the simulated data,
# with varying levels of measurement noise
#

source("common.R")

# import simulated MIDs, time point 75
sim_mi_data <- readRDS(file.path(mi_data_path, 'simulated_mi_data_t75.rds'))
n_experiments <- length(sim_mi_data$experiments)

# add dirichlet noise and compute distance matrices

stdevs <- c(0.01, 0.05, 0.1)
n_resamples <- 10

{
   set.seed(7462931)
   
   for (stdev in stdevs) {
      create_dir_if_not_exists(sim_dir_path(stdev))
      for (rep_nr in 1:n_resamples) {
         print("add noise")
         sim_mi_data_noisy <- midata_transform(sim_mi_data, function(x)
            c(random_mid(x, stdev, 1)))
         saveRDS(sim_mi_data_noisy, sim_mi_data_path(stdev, rep_nr))
         
         print("compute dm, all experiments")
         assign_list[sim_dm, sim_conv_index] <- conv_reduce_all(
            sim_mi_data_noisy,
            1:n_experiments,
            f = midist::euclidean_sum_dist,
            g = which.min
         )
         sim_dm <- replace_na_with_max(sim_dm)
         saveRDS(sim_dm, sim_dm_path(stdev, rep_nr))
         
         if(stdev == 0.01) {
            print("compute dm, individual experiments")
            for(experiment in sim_mi_data$experiments) {
               cat(experiment, "")
               assign_list[sim_dm, sim_conv_index] <- conv_reduce_all(
                  sim_mi_data_noisy,
                  experiment,
                  f = midist::euclidean_dist,
                  g = which.min
               )
               sim_dm <- replace_na_with_max(sim_dm)
               saveRDS(sim_dm, sim_dm_ind_exp_path(stdev, experiment, rep_nr))
            }
         }
      }
   }
}

# without noise

{
   create_dir_if_not_exists(sim_dir_path(stdev = 0))
   assign_list[sim_dm, sim_conv_index] <- conv_reduce_all(
      sim_mi_data,
      1:n_experiments,
      f = midist::euclidean_sum_dist,
      g = which.min
   )
   sim_dm <- replace_na_with_max(sim_dm)
   saveRDS(sim_dm, sim_dm_path(stdev = 0, rep_nr = 1))
}

#
# distance matrices for subsets of metabolites
#

stdev <- 0.01
subset_sizes <- c(300, 200, 100, 50)
n_resamples <- 10

{
   set.seed(4725116)
   
   for (subset_size in subset_sizes) {
      create_dir_if_not_exists(subset_dir_path(subset_size))
      for (rep_nr in 1:n_resamples) {
         print("subsample and add noise")
         sim_mi_data_subset <- sim_mi_data %>%
            midata_subset(sample(sim_mi_data$peak_ids, subset_size)) %>%
            midata_transform(function(x) c(random_mid(x, stdev, 1)))
         saveRDS(sim_mi_data_subset, subset_mi_data_path(subset_size, rep_nr))
         
         print("compute dm")
         assign_list[sim_dm, sim_conv_index] <- conv_reduce_all(
            sim_mi_data_subset,
            1:n_experiments,
            f = midist::euclidean_sum_dist,
            g = which.min
         )
         sim_dm <- replace_na_with_max(sim_dm)
         saveRDS(sim_dm, subset_dm_path(subset_size, rep_nr))
      }
   }
}

