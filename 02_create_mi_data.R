#
#  Create MI data object
#

library(remn)

source("common.R")

# Read HMEC peak areas
hmec_peak_list <- read_hmec_peak_list()
hmec_peak_areas <- read_hmec_peak_areas()

experiment_names <- substr(
    colnames(hmec_peak_areas)[3:length(hmec_peak_areas)],
    1, 3
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
