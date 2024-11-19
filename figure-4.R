#
# Figure 4: Systematic peak annotation & validation
#

source("common.R")


hmec_peak_list <- read_hmec_peak_list()
peak_hmdb_compound <- read_peak_hmdb_compound()
hmec_mi_data <- readRDS(file.path(mi_data_path, "hmec_mi_data_censored.rds"))
hmec_dm <- readRDS(file.path(mid_distance_path, 'hmec_dm.rds'))
plotly_tooltips <- read_plotly_tooltips()


#
# Figure 4a statistics
#



# number of unknown peaks, lipids and non-lipids
hmec_peak_list %>%
    mutate(is_known = is.na(known_met_id)) %>%
    group_by(is_known, is_lipid) %>% count()

#
# Figure 4b, 4889 aspartyl-glucosamine
#

example_peak_id <- "4889"
n_neighbors <- 20

# candidate annotations for this peak
peak_hmdb_compound %>% filter(peak_id == example_peak_id)

# nearest neighbors by MID distance
neighbors <- names(sort(hmec_dm[example_peak_id, ])[1:n_neighbors])

umap_proj <- umap_projection(
    hmec_dm[neighbors, neighbors], n_neighbors = 5, random_seed = 35179311)

plot_umap(umap_proj)

plot_umap_interactive(umap_proj, plotly_tooltips[neighbors,"tooltip"])

# list of known metabolite ids in this neighborhood
nearest_known_table(neighbors, hmec_peak_list)


#
# Figure 4c number of predictions
#

# no predicted structure for lipids
hmec_peak_list %>% filter((is_lipid == 1) & !is.na(inchi_key)) %>% nrow

# number of new annotations (excluding previously known)
hmec_peak_list %>%
    filter((is_lipid == 0) & is.na(known_met_id) & !is.na(inchi_key)) %>% nrow
# of which new unique strucures
hmec_peak_list %>%
    filter((is_lipid == 0) & is.na(known_met_id) & !is.na(inchi_key)) %>% 
    distinct(inchi_key) %>% nrow

# GNPS annotations
gnps_annotations <- read_tsv(file.path(input_data_path, "gnps-annotations.tsv")) %>%
    mutate(peak_id = as.character(peak_id)) %>%
    mutate(skeleton = str_sub(inchi_key, end = 14))

# verify peak ids are valid
stopifnot(
    gnps_annotations %>%
        anti_join(hmec_peak_list, join_by('peak_id')) %>% nrow() == 0
)

# number of GNPS-annotated peaks, lipids and non-lipids
gnps_annotations %>%
    inner_join(
        hmec_peak_list %>% select(peak_id, is_lipid),
        join_by('peak_id')) %>%
    group_by(is_lipid) %>% count()

# number of unique compounds
gnps_annotations %>%
    inner_join(
        hmec_peak_list %>% select(peak_id, is_lipid),
        join_by('peak_id')) %>%
    distinct(is_lipid, skeleton) %>%
    group_by(is_lipid) %>% count()


#
# Figure 4d spermidine-related compounds
#


# neighborhood
top_n <- 10
top_index <- order(hmec_dm["6174", ])[1:top_n]

umap_proj <- umap_projection(
    hmec_dm[top_index, top_index], n_neighbors = 5, random_seed = 759284725)
plot_umap(umap_proj)
plot_umap_interactive(umap_proj, plotly_tooltips$tooltip[top_index])

# MID plots
selected_exp <- c("glc", "gln", "arg", "met")

# 6174 C10 compound
plot_mid_matrix(
    c13correct_cols(
        get_mid_matrix(hmec_mi_data, "6174", selected_exp)
    ),
    max_mi_fraction = 0.3
)

# 6284 C9 compound
plot_mid_matrix(
   c13correct_cols(
      get_mid_matrix(hmec_mi_data, "6284", selected_exp)
   ),
   max_mi_fraction = 0.3
)


# write_mid_matrix_image(
#     "C:/tmp/6335-polyamine.png",
#     get_mid_matrix(hmec_mi_data, "6335", selected_exp),
#     max_mi_fraction = 0.3)


#
# Fig 4e 13c enrichmennt for glutamine co-eluting peaks
#

selected_exp <- c("arg", "asn", "gln")
coeluting_peaks <- rev(c("4905", "5991", "6698", "6849", "7072"))

# 13C enrichment for each peak and experiment
enrichment_mat <- sapply(
    coeluting_peaks,
    function(p) {
       apply(
           get_avg_mid(hmec_mi_data, p, selected_exp), 2,
           isotopic_enrichment
       )
    }
)

image(
    enrichment_mat,
    col = colorRampPalette(c("white", "blue"))(100)
)


