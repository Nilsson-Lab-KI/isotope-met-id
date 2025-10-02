
source("common.R")


hmec_mi_data <- readRDS(file.path(mi_data_path, 'hmec_mi_data_censored.rds'))
n_peaks <- length(hmec_mi_data$peak_ids)
n_experiments <- length(hmec_mi_data$experiments)

# compute distance matrix using all experiments except unlabeled
# we impute size 1 metabolites only (notably CO2)
# all other missing convolutions will yield NA distance
assign_list[hmec_dm, hmec_conv_index] <- conv_reduce_all(
    hmec_mi_data, 1:(n_experiments - 1),
    f = midist::euclidean_sum_dist,
    g = which.min,
    impute = 1
)

# impute missing values with maximal distance
max_distance = max(hmec_dm, na.rm = TRUE)
hmec_dm[which(is.na(hmec_dm))] <- max_distance

# Write distance matrix
saveRDS(
    hmec_dm,
    file.path(mid_distance_path, 'hmec_dm.rds')
)
saveRDS(
    hmec_conv_index,
    file.path(mid_distance_path, 'hmec_conv_index.rds')
)


# same, for 13C-corrected data
hmec_mi_data_13c_corr <- readRDS(file.path(mi_data_path, 'hmec_mi_data_13c_corr.rds'))

assign_list[hmec_13c_corr_dm, hmec_13c_corr_conv_index] <- conv_reduce_all(
   hmec_mi_data_13c_corr, 1:(n_experiments - 1),
   f = midist::euclidean_sum_dist,
   g = which.min,
   impute = 1
)

# impute missing values with maximal distance
max_13c_corr_distance = max(hmec_13c_corr_dm, na.rm = TRUE)
hmec_13c_corr_dm[which(is.na(hmec_13c_corr_dm))] <- max_13c_corr_distance

# Write distance matrix
saveRDS(
   hmec_13c_corr_dm,
   file.path(mid_distance_path, 'hmec_13c_corr_dm.rds')
)
saveRDS(
   hmec_13c_corr_conv_index,
   file.path(mid_distance_path, 'hmec_13c_corr_conv_index.rds')
)

