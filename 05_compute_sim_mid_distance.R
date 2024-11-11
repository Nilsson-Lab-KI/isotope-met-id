#
# Compute MID distances matrices from the simulated data,
# with varying levels of measurement noise
#

source("common.R")

# TODO: this rds file was computed previousl
# need to regenerate it in 02_create_mi_data.R
sim_mi_data <- readRDS(file.path(mi_data_path, 'sim_mi_data.rds'))
n_experiments <- length(sim_mi_data$experiments)

# add dirichlet noise, compute distance matrices and repeat

stdevs <- c(0.01, 0.05, 0.1)
n_resamples <- 1

{
   set.seed(7462931)
   
   for (stdev in stdevs) {
      create_dir_if_not_exists(sim_dir_path(stdev))
      for (rep_nr in 1:n_resamples) {
         print("add noise")
         sim_mi_data_noisy <- midata_transform(sim_mi_data, function(x)
            c(random_mid(x, stdev, 1)))
         saveRDS(sim_mi_data_noisy, sim_mi_data_path(stdev, rep_nr))
         
         print("compute dm")
         assign_list[sim_dm, sim_conv_index] <- conv_reduce_all(
            sim_mi_data_noisy,
            1:n_experiments,
            f = midist::euclidean_sum_dist,
            g = which.min
         )
         # impute missing values with maximal distance
         max_distance <- max(sim_dm, na.rm = TRUE)
         sim_dm[which(is.na(sim_dm))] <- max_distance
         
         saveRDS(sim_dm, sim_dm_path(stdev, rep_nr))
      }
   }
}
