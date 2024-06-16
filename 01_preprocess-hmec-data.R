#
# Preprocess input data
#

source("common.R")


# preprocess HMEC data

hmec_peak_list <- read_hmec_peak_list()

# read peak areas for all peaks
peak_areas_path <- file.path(input_data_path, "hmec-peak-areas.tsv")
hmec_peak_areas <- read_tsv(peak_areas_path)

# order tracer data columns alphabetic by tracer, unlabeled data last
tracer_data_col_names = colnames(hmec_peak_areas)[3:65]
column_order <- c(1:2, 2 + order(tracer_data_col_names), 66:68)

# subset to QC peaks
hmec_peak_areas <- hmec_peak_areas[
    which(hmec_peak_areas$peak_id %in% hmec_peak_list$peak_id),
    column_order
]

# write QC peak areas
write.table(
    hmec_peak_areas,
    file.path(preprocessed_data_path, "qc-peak-areas.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
)


# create tooltip strings for interactive plots

peak_hmdb_compound <- read_peak_hmdb_compound()

plotly_tooltips <- data.frame(peak_id = hmec_peak_list$peak_id) %>%
    left_join(
        peak_hmdb_compound,
        by = 'peak_id') %>%
    # type conversion to get the correct sort order
    mutate(peak_id = as.integer(peak_id)) %>%
    group_by(peak_id) %>%
    summarize(tooltip = str_c(name, collapse = "|")) %>%
    mutate(peak_id = as.character(peak_id))

stopifnot(
    all(plotly_tooltips$peak_id == hmec_peak_list$peak_id)
)

write.table(
    plotly_tooltips,
    file.path(preprocessed_data_path, "plotly_tooltips.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
)


