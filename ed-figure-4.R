#
#  Extended Data Fig 4
#

source("common.R")

hmec_peak_list <- read_hmec_peak_list()
peak_hmdb_compound <- read_peak_hmdb_compound()
hmec_mi_data <- readRDS(file.path(mi_data_path, "hmec_mi_data_censored.rds"))
hmec_dm <- readRDS(file.path(mid_distance_path, 'hmec_dm.rds'))
plotly_tooltips <- read_plotly_tooltips()


#
# ED Figure 4a lipid MIDs
#

max_mid_length <- max(hmec_mi_data$peak_n_atoms) + 1

pad_right <- function(x, n)
{
    c(x, rep(0, n - length(x)))
}

lipid_n_carbons <- hmec_mi_data$peak_n_atoms[which(hmec_peak_list$is_lipid == 1)]

lipid_glc_mids <- apply(
    hmec_peak_list %>% filter(is_lipid == 1) %>% select(peak_id),
    MARGIN = 1,
    function(peak_id) pad_right(
        get_avg_mid(hmec_mi_data, peak_id, "glc"),
        max_mid_length
    )
)

plot_mid_matrix(
    c13correct_cols(lipid_glc_mids[, order(lipid_n_carbons)]),
    max_mi_fraction = 0.3
)

plot(lipid_n_carbons[order(lipid_n_carbons)], type='l')

# write_mid_matrix_image(
#     'C:/tmp/lipid_mids.png',
#     lipid_glc_mids[, order(lipid_n_carbons)],
#     max_mi_fraction = 0.3
# )


#
# ED Figure 4b
# number of candidate structures per peak, for polar compounds
#

n_candidates = peak_hmdb_compound %>%
    inner_join(
        hmec_peak_list %>% filter(is_lipid == 0) %>% select(peak_id),
        by = "peak_id"
    ) %>%
    group_by(peak_id) %>% count()

n_candidates$n %>% table %>% barplot

# average
mean(n_candidates$n)


#
# ED Figure 4c, 5579 nicotinamide riboside
#

example_peak_id <- "5579"
n_neighbors <- 20

# candidate annotations for this peak
peak_hmdb_compound %>% filter(peak_id == example_peak_id)

# nearest neighbors by MID distance
neighbors <- names(sort(hmec_dm[example_peak_id, ])[1:n_neighbors])

umap_proj <- umap_projection(
    hmec_dm[neighbors, neighbors], n_neighbors = 5, random_seed = 82735681)

plot_umap(umap_proj)

plot_umap_interactive(umap_proj, plotly_tooltips[neighbors,"tooltip"])

# list of known metabolite ids in this neighborhood
nearest_known_table(neighbors, hmec_peak_list)


#
# ED Figure 4d, 2147 N-acetyl-threonine
#

example_peak_id <- "2147"
n_neighbors <- 20

# nearest neighbors by MID distance
neighbors <- names(sort(hmec_dm[example_peak_id, ])[1:n_neighbors])

umap_proj <- umap_projection(
    hmec_dm[neighbors, neighbors], n_neighbors = 5, random_seed = 4132512)

plot_umap(umap_proj)

plot_umap_interactive(umap_proj, plotly_tooltips[neighbors,"tooltip"])

# list of known metabolite ids in this neighborhood
nearest_known_table(neighbors, hmec_peak_list)


#
# ED Figure 4e, 5128 NAAG
#

example_peak_id <- "5128"
n_neighbors <- 20

# nearest neighbors by MID distance
neighbors <- names(sort(hmec_dm[example_peak_id, ])[1:n_neighbors])

umap_proj <- umap_projection(
    hmec_dm[neighbors, neighbors], n_neighbors = 5, random_seed = 82735681)

plot_umap(umap_proj)

plot_umap_interactive(umap_proj, plotly_tooltips[neighbors,"tooltip"])

# list of known metabolite ids in this neighborhood
nearest_known_table(neighbors, hmec_peak_list)

plot_mid_matrix(
    c13correct_cols(
        get_mid_matrix(hmec_mi_data, example_peak_id, c("asp", "glc", "gln", "glu"))
    ),
    max_mi_fraction = 0.3
)

#
# ED Figure 4j, 5872 unknown histidine compound
#
example_peak_id <- "5872"
n_neighbors <- 20

# none of the HMDB annotations are related to histidine
peak_hmdb_compound %>% filter(peak_id == example_peak_id)

# nearest neighbors by MID distance
neighbors <- names(sort(hmec_dm[example_peak_id, ])[1:n_neighbors])

umap_proj <- umap_projection(
    hmec_dm[neighbors, neighbors], n_neighbors = 5, random_seed = 759284725)

plot_umap(umap_proj)

plot_umap_interactive(umap_proj, plotly_tooltips[neighbors,"tooltip"])

# list of known metabolite ids in this neighborhood
nearest_known_table(neighbors, hmec_peak_list)


#
# ED Figure 4, 5450 glucosyl-serine
#

example_peak_id <- "5450"
n_neighbors <- 20

peak_hmdb_compound %>% filter(peak_id == example_peak_id)

# nearest neighbors by MID distance
neighbors <- names(sort(hmec_dm[example_peak_id, ])[1:n_neighbors])

umap_proj <- umap_projection(
    hmec_dm[neighbors, neighbors], n_neighbors = 5, random_seed = 23218353)

plot_umap(umap_proj)

plot_umap_interactive(umap_proj, plotly_tooltips[neighbors,"tooltip"])

# list of known metabolite ids in this neighborhood
nearest_known_table(neighbors, hmec_peak_list)

plot_mid_matrix(
    c13correct_cols(
        get_mid_matrix(hmec_mi_data, example_peak_id, c("ser", "glc"))
    ),
    max_mi_fraction = 0.5
)

# write_mid_matrix_image(
#     "C:/tmp/5450-glucosyl-serine-mid.png",
#     get_mid_matrix(hmec_mi_data, "5450", c("ser", "glc")),
#     max_mi_fraction = 0.5)


#
# ED Figure 4m, 6133 lactoyl-asparagine
#
example_peak_id <- "6133"
n_neighbors <- 20

peak_hmdb_compound %>% filter(peak_id == example_peak_id)

# nearest neighbors by MID distance
neighbors <- names(sort(hmec_dm[example_peak_id, ])[1:n_neighbors])

umap_proj <- umap_projection(
    hmec_dm[neighbors, neighbors], n_neighbors = 5, random_seed = 96529411)

plot_umap(umap_proj)

plot_umap_interactive(umap_proj, plotly_tooltips[neighbors,"tooltip"])

# list of known metabolite ids in this neighborhood
nearest_known_table(neighbors, hmec_peak_list)

plot_mid_matrix(
    c13correct_cols(
        get_mid_matrix(hmec_mi_data, example_peak_id, c("asn", "glc"))
    ),
    max_mi_fraction = 0.5
)

# write_mid_matrix_image(
#     "C:/tmp/6133-lactoyl-asparagine-mid.png",
#     get_mid_matrix(hmec_mi_data, "6133", c("asn", "glc")),
#     max_mi_fraction = 0.5)

