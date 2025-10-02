#
#  Create MI data objects
#

source("common.R")

#
# HMEC data
#

hmec_peak_areas <- read_hmec_peak_areas()

# drop the experiment replicate number from column names
experiment_names <- colnames(hmec_peak_areas)[3:length(hmec_peak_areas)]
experiment_names <- substr(
    experiment_names,
    1, nchar(experiment_names)-2
)


# Create MIData object
hmec_mi_data <- MIData(
    hmec_peak_areas,
    exp_names = experiment_names,
    exp_columns = 3:length(hmec_peak_areas)
)

n_peaks <- length(hmec_mi_data$peak_ids)
n_experients <- length(hmec_mi_data$experiments)

# Write MIData object
saveRDS(
    hmec_mi_data,
    file.path(mi_data_path, 'hmec_mi_data.rds')
)

# create corresponding text file for Suppl Data 2
get_mid_df <- function(mi_data, peak_id)
{
    mid_df <- as.data.frame(
        get_mid_matrix(mi_data, peak_id, mi_data$experiments))
    # TODO: this should go into the remn package
    colnames(mid_df) <- unlist(
        lapply(
            1:length(mi_data$experiments),
            function(i)
                paste(
                    rep(mi_data$experiments[i], mi_data$exp_n_rep[i]),
                    1:mi_data$exp_n_rep[i],
                    sep = '_'
                )
        )
    )
    return(
        mid_df %>%
            mutate(peak_id = peak_id) %>% relocate(peak_id) %>%
            mutate(mi = 0:(nrow(mid_df)-1)) %>% relocate(mi, .after = peak_id)
    )
}

hmec_all_mids <- lapply(
    hmec_mi_data$peak_ids,
    function(peak_id) get_mid_df(hmec_mi_data, peak_id)) %>% bind_rows()

write.table(
    hmec_all_mids, file.path(mi_data_path, 'suppl_data_2.tsv'),
    sep = '\t', row.names = FALSE
)

# export a table of carbon numbers
hmec_peak_n_carbons <- data.frame(
   peak_id = hmec_mi_data$peak_ids,
   n_carbons = hmec_mi_data$peak_n_atoms
)

write.table(
   hmec_peak_n_carbons,
   file.path(preprocessed_data_path, "hmec_peak_n_carbons_qc.tsv"),
   sep = "\t", row.names = FALSE, quote = FALSE
)


# Censor false MIs
hmec_mi_data_censored <- censor_false_mi(
    hmec_mi_data,
    min_experiments = 10
)

# Write MIData object
saveRDS(
    hmec_mi_data_censored,
    file.path(mi_data_path, 'hmec_mi_data_censored.rds')
)


#
# Natural 13C correction
#

hmec_mi_data_13c_corr <- correct_natural_13c(hmec_mi_data_censored)

saveRDS(
   hmec_mi_data_13c_corr,
   file.path(mi_data_path, 'hmec_mi_data_13c_corr.rds')
)


#
# Simulated data
#

sim_peak_areas <- read_tsv(
   file.path(input_data_path, "simulated_mids_t75.tsv"))

experiment_names <- colnames(sim_peak_areas)[-(1:2)]

# replicate each experiment 3 times and create a new data frame
sim_peak_areas_replicated <- cbind(
   sim_peak_areas[1:2],
      do.call(
      cbind,
      lapply(
         experiment_names,
         function(exp_name) {
            sapply(
               paste(exp_name, 1:3, sep="_"),
               function(x) sim_peak_areas[[exp_name]])
         }
      )
   )
)

experiment_names_replicated <- sim_peak_areas_replicated[-(1:2)] %>%
   colnames %>% str_sub(end = -3)

# Create MIData object
sim_mi_data <- MIData(
   sim_peak_areas_replicated,
   exp_names = experiment_names_replicated,
   exp_columns = 3:length(sim_peak_areas_replicated)
)

# Write MIData object
saveRDS(
   sim_mi_data,
   file.path(mi_data_path, 'simulated_mi_data_t75.rds')
)

