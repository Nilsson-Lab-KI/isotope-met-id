#
# Figure 5
#

source("common.R")


hmec_peak_list <- read_hmec_peak_list()
peak_hmdb_compound <- read_peak_hmdb_compound()
hmec_mi_data <- readRDS(file.path(mi_data_path, "hmec_mi_data_censored.rds"))
hmec_dm <- readRDS(file.path(mid_distance_path, 'hmec_dm.rds'))
plotly_tooltips <- read_plotly_tooltips()


#
# Fig 5a, peak 5665 (TMGL) neighborhood local projection
#

example_peak_id <- "5665"
n_neighbors <- 30

# candidate annotations for this peak
peak_hmdb_compound %>% filter(peak_id == example_peak_id)

# nearest neighbors by MID distance
neighbors <- names(sort(hmec_dm[example_peak_id, ])[1:n_neighbors])

umap_proj <- umap_projection(
   hmec_dm[neighbors, neighbors], n_neighbors = 7, random_seed = 35261820)

plot_umap(umap_proj)

plot_umap_interactive(umap_proj, plotly_tooltips[neighbors, "tooltip"])

# list of known metabolite ids in this neighborhood
nearest_known_table(neighbors, hmec_peak_list)


#
# Fig 5b TMGL and related MIDs
#
selected_exp <- c("gly", "lys", "met", "ser")

plot_mid_matrix(
   c13correct_cols(
      get_mid_matrix(hmec_mi_data, "5665", selected_exp)
   ),
   max_mi_fraction = 0.3
)

# 6269 trimethyllysine +H (known)
plot_mid_matrix(
   c13correct_cols(
      get_mid_matrix(hmec_mi_data, "6269", selected_exp)
   ),
   max_mi_fraction = 0.3
)

# 2531 glycine (known)
plot_mid_matrix(
   c13correct_cols(
      get_mid_matrix(hmec_mi_data, "2531", selected_exp)
   ),
   max_mi_fraction = 0.3
)

# convolution MID tmlys * gly
mids_tmlys_gly_conv <- convolute_all(
   get_mid_matrix(hmec_mi_data, "6269", selected_exp),
   get_mid_matrix(hmec_mi_data, "2531", selected_exp)
)
plot_mid_matrix(
   c13correct_cols(mids_tmlys_gly_conv),
   max_mi_fraction = 0.3
)

