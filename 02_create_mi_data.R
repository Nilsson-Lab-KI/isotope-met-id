#
#  Create MI data object
#

library(remn)

source("common.R")

# Read HMEC peak areas
hmec_peak_list <- read_hmec_peak_list()
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
