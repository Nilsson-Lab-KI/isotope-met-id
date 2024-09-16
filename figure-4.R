#
# Figure 4: Systematic peak annotation & validation
#

source("common.R")

library(ggplot2)
library(plotly)
library(remn)
library(stringr)


hmec_peak_list <- read_hmec_peak_list()

peak_hmdb_compound <- read_peak_hmdb_compound()

hmec_mi_data <- readRDS(file.path(mi_data_path, "hmec_mi_data_censored.rds"))

hmec_dm <- readRDS(file.path(mid_distance_path, 'hmec_dm.rds'))

plotly_tooltips <- read_plotly_tooltips()


#
# Figure 4b, 1583 glycerol-P-glycerol
#

# nearest neighbors by MID distance
top_index <- order(hmec_dm["1583", ])[1:20]

umap_proj <- umap_projection(
    hmec_dm[top_index, top_index], n_neighbors = 5, random_seed = 3614793)

plot_umap(umap_proj)

plot_umap_interactive(umap_proj, plotly_tooltips$tooltip[top_index])


#
# ED Figure 4a
# 

# number of candidate structures per peak, for polar compounds only
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
# ED Figure 4b lipid MIDs
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

write_mid_matrix_image(
    'C:/tmp/lipid_mids.png',
    lipid_glc_mids[, order(lipid_n_carbons)],
    max_mi_fraction = 0.3
)


#
# Figure 4c, 4889 aspartyl-glucosamine
#

peak_hmdb_compound %>% filter(peak_id == 4889)

# nearest neighbors by MID distance
top_index <- order(hmec_dm["4889", ])[1:20]

umap_proj <- umap_projection(
    hmec_dm[top_index, top_index], n_neighbors = 5, random_seed = 35179311)

plot_umap(umap_proj)
plot_umap_interactive(umap_proj, plotly_tooltips$tooltip[top_index])

plot_mid_matrix(
    c13correct_cols(
        get_mid_matrix(hmec_mi_data, "4889", c("asn", "glc"))
    ),
    max_mi_fraction = 0.3
)

write_mid_matrix_image(
    "C:/tmp/aspartyl-glucosamine-mid.png",
    get_mid_matrix(hmec_mi_data, "4889", c("asn", "glc")),
    max_mi_fraction = 0.3)


#
# Figure 4e 6174 spermic acid? 6284 acetyl-spermidine
#

# nearest neighbors by MID distance
top_index <- order(hmec_dm["6174", ])[1:10]

umap_proj <- umap_projection(
    hmec_dm[top_index, top_index], n_neighbors = 5, random_seed = 759284725)

plot_umap(umap_proj)
plot_umap_interactive(umap_proj, plotly_tooltips$tooltip[top_index])

# candidates
peak_hmdb_compound %>% filter(peak_id == 6174)

# MID plots
selected_exp <- c("glc", "gln", "arg", "met")

plot_mid_matrix(
    c13correct_cols(
        get_mid_matrix(hmec_mi_data, "6174", selected_exp)
    ),
    max_mi_fraction = 0.3
)

write_mid_matrix_image(
    "C:/tmp/6174-polyamine.png",
    get_mid_matrix(hmec_mi_data, "6174", selected_exp),
    max_mi_fraction = 0.3)


plot_mid_matrix(
    c13correct_cols(
        get_mid_matrix(hmec_mi_data, "6335", selected_exp)
    ),
    max_mi_fraction = 0.1
)


#
# ED Figure 4i, 6174 unknown histidine compound
#
peak_hmdb_compound %>% filter(peak_id == 5872)

plot_mid_matrix(
    c13correct_cols(
        get_mid_matrix(hmec_mi_data, "5872", hmec_mi_data$experiments)
    ),
    max_mi_fraction = 0.3
)

# nearest neighbors by MID distance
top_index <- order(hmec_dm["5872", ])[1:20]

umap_proj <- umap_projection(
    hmec_dm[top_index, top_index], n_neighbors = 5, random_seed = 759284725)

plot_umap(umap_proj)
plot_umap_interactive(umap_proj, plotly_tooltips$tooltip[top_index])


#
# ED Figure 4j, 5450 glucosyl-serine
#
peak_hmdb_compound %>% filter(peak_id == 5450)

plot_mid_matrix(
    c13correct_cols(
        get_mid_matrix(hmec_mi_data, "5450", c("ser", "glc"))
    ),
    max_mi_fraction = 0.5
    
)

write_mid_matrix_image(
    "C:/tmp/5450-glucosyl-serine-mid.png",
    get_mid_matrix(hmec_mi_data, "5450", c("ser", "glc")),
    max_mi_fraction = 0.5)

# nearest neighbors by MID distance
top_index <- order(hmec_dm["5450", ])[1:20]

umap_proj <- umap_projection(
    hmec_dm[top_index, top_index], n_neighbors = 5, random_seed = 23218353)

plot_umap(umap_proj)
plot_umap_interactive(umap_proj, plotly_tooltips$tooltip[top_index])


#
# ED Figure 4k, 61330 lactoyl-asparagine
#
peak_hmdb_compound %>% filter(peak_id == 6133)

plot_mid_matrix(
    c13correct_cols(
        get_mid_matrix(hmec_mi_data, "6133", c("asn", "glc"))
    ),
    max_mi_fraction = 0.5
    
)

write_mid_matrix_image(
    "C:/tmp/6133-lactoyl-asparagine-mid.png",
    get_mid_matrix(hmec_mi_data, "6133", c("asn", "glc")),
    max_mi_fraction = 0.5)

# nearest neighbors by MID distance
top_index <- order(hmec_dm["6133", ])[1:20]

umap_proj <- umap_projection(
    hmec_dm[top_index, top_index], n_neighbors = 5, random_seed = 96529411)

plot_umap(umap_proj)
plot_umap_interactive(umap_proj, plotly_tooltips$tooltip[top_index])


#
# ED Figur 4x, 6505 spinacine
#
peak_hmdb_compound %>% filter(peak_id == 6505)

plot_mid_matrix(
    c13correct_cols(
        get_mid_matrix(hmec_mi_data, "6505", c("his"))
    ),
    max_mi_fraction = 0.3
)

write_mid_matrix_image(
    "C:/tmp/spinacine-mid.png",
    get_mid_matrix(hmec_mi_data, "6505", c("his")),
    max_mi_fraction = 1.0)


#
# ED Figure 4c, 5579 nicotinamide riboside
#

peak_hmdb_compound %>% filter(peak_id == 5579)

# nearest neighbors by MID distance
top_index <- order(hmec_dm["5579", ])[1:20]

which(names(hmec_dm["5579", ]) == "1679")

umap_proj <- umap_projection(
    hmec_dm[top_index, top_index], n_neighbors = 5, random_seed = 82735681)

plot_umap(umap_proj)
plot_umap_interactive(umap_proj, plotly_tooltips$tooltip[top_index])

plot_mid_matrix(
    c13correct_cols(
        get_mid_matrix(hmec_mi_data, "5579", c("g"))
    ),
    max_mi_fraction = 0.3
)

#
# ED Figure 4d, 2146 N-acetyl-threonine
#

peak_hmdb_compound %>% filter(peak_id == 2147)

# nearest neighbors by MID distance
top_index <- order(hmec_dm["2147", ])[1:20]

umap_proj <- umap_projection(
    hmec_dm[top_index, top_index], n_neighbors = 5, random_seed = 4132512)

plot_umap(umap_proj)
plot_umap_interactive(umap_proj, plotly_tooltips$tooltip[top_index])


#
# ED Figure 4e, 5128 NAAG
#

peak_hmdb_compound %>% filter(peak_id == 5128)

# nearest neighbors by MID distance
top_index <- order(hmec_dm["5128", ])[1:20]

umap_proj <- umap_projection(
    hmec_dm[top_index, top_index], n_neighbors = 5, random_seed = 82735681)

plot_umap(umap_proj)
plot_umap_interactive(umap_proj, plotly_tooltips$tooltip[top_index])

zplot_mid_matrix(
    c13correct_cols(
        get_mid_matrix(hmec_mi_data, "5128", c("asp", "glc", "gln", "glu"))
    ),
    max_mi_fraction = 0.3
)


# examples

# 6710 burytobetaine
peak_hmdb_compound %>% filter(peak_id == 6710)
# 1257 NAAG
peak_hmdb_compound %>% filter(peak_id == 1257)
# 1583 glycerol-P-glycerol
peak_hmdb_compound %>% filter(peak_id == 1583)
# 5667 isovaleryl-carnitine
peak_hmdb_compound %>% filter(peak_id == 5667)
# 5656 lactoyl-arginine
peak_hmdb_compound %>% filter(peak_id == 5656)



#
# nearest neighbor based peak annotation
#

# Create MIData object
hmec_mi_data <- readRDS(file.path(mi_data_path, 'hmec_mi_data_censored.rds'))

n_peaks <- length(hmec_mi_data$peak_ids)
n_experiments <- length(hmec_mi_data$experiments)

# verify that MIData peak IDs matches peak list
stopifnot(all(hmec_mi_data$peak_ids == hmec_peak_list$peak_id))

# distance matrix based on all experiments except unlabeled
assign_list[dm_with_na, conv_index] <- conv_reduce_all(
    hmec_mi_data, 1:(n_experiments - 1),
    f = remn::euclidean_sum_dist,
    g = which.min,
    impute = 1
)
# replace NA values with largest distance
dm <- dm_with_na
dm[which(is.na(dm))] <- max(dm_with_na, na.rm = TRUE)
hist(dm[lower.tri(dm)], n = 100)

#
# Nearest known for each unknown
#

# known peaks

# known_peak_ids <- peak_annotations %>%
#     filter(known_met_id != "") %>% { .$peak_id } %>% as.character
# 
# find_nearest <- function(dm, peak, known_peaks, n)
# {
#     top_known <- sort(dm[peak, known_peaks])[2:(n+1)]
#     return(
#         data.frame(
#             peak_id = names(top_known),
#             mid_dist = top_known
#         )
#     )
# }
# 
# find_nearest(dm, "1583", known_peak_ids, 10) %>%
#     inner_join(
#         peak_annotations %>% select(peak_id, known_met_id),
#         by = "peak_id"
#     )


#unit_vector <- function(n, k) ifelse(1:n == k, 1, 0)

# find_nearest_known <- function(dm, conv_index, known_ids)
# {
#     known_index <- match(known_ids, rownames(dm))
#     unknown_index <- setdiff(1:nrow(dm), known_index)
#     best_known <- sapply(
#         unknown_index,
#         function(k) known_index[which.min(dm[k, known_index])]
#     )
#     best_known_conv <- conv_index[cbind(unknown_index, best_known)]
#     data.frame(
#         unknown_id = as.integer(rownames(dm)[unknown_index]),
#         known_id = as.integer(rownames(dm)[best_known]),
#         conv_id = as.integer(rownames(dm)[best_known_conv]),
#         distance = dm[cbind(unknown_index, best_known)]
#     )
# }
# 
# known_ids <- hmec_peak_list %>%
#     filter(is_known == 1) %>% { as.character(.$peak_id) }
# 
# best_known <- find_nearest_known(dm, conv_index, known_ids)
# 
# best_known <- best_known %>%
#     inner_join(
#         peak_annotations %>% select(peak_id, hypothesis),
#         by = join_by(unknown_id == peak_id)
#     ) %>%
#     inner_join(
#         peak_annotations %>% select(peak_id, known_met_id),
#         by = join_by(known_id == peak_id)
#     ) %>%
#     left_join(
#         peak_annotations %>% select(peak_id, hypothesis) %>%
#             rename(conv_hypothesis = hypothesis),
#         by = join_by(conv_id == peak_id)
#     )


#
# 6710 butyro-betaine
#

# selected experiments to plot
selected_exp <- c("lys", "met")

mids_6710 <- get_mid_matrix(hmec_mi_data, "6710", selected_exp)
plot_mid_matrix(c13correct_cols(mids_6710)[-1, ])


#
# Glutamine co-eluting peaks
#

selected_exp <- c("arg", "asn", "gln")
coeluting_peaks <- c("7072", "6849", "6698", "5991", "4905")

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

write_image(
    file.path(data_path, "plots", "HMEC-enrichment-clustered-heatmap.png"),
    reverse_rows(t(enrichment_mat[enrichment_clust$order, ])),
    colorRamp(c("white", "blue"))
)


#
# 5665 TMLG neighborhood local projection
#

# HMDB candidates
peak_skeleton %>% filter(peak_id == 5665)

# nearest known
find_nearest(dm, "5665", known_peak_ids, 10) %>%
    inner_join(
        peak_annotations %>% select(peak_id, known_met_id),
        by = "peak_id"
    )
# convolutant for TML
hmec_peak_list[conv_index["6269", "5665"], ]
 

# nearest neighbors
top_index <- order(dm["5665", ])[1:30]
top_neighbors <- hmec_peak_list[top_index, ]

# fix random seed for reproducible UMap projection
set.seed(35261820)

umap_config <- umap.defaults
umap_config$n_neighbors <- 7
umap_proj <- umap(dm[top_index, top_index], input = "dist", config = umap_config)

umap_df <- data.frame(
    peak_id = hmec_peak_list$peak_id[top_index],
    peak_ann = paste(hmec_peak_list$peak_id, hmec_peak_list$netid_annotation)[top_index],
    umap_1 = umap_proj$layout[, 1],
    umap_2 = umap_proj$layout[, 2]
)

# plot projection
umap_df %>%
    ggplot(aes(x = umap_1, y = umap_2)) +
    geom_point(alpha = 0.7, colour = "black") +
    theme_classic()

# interactive plot with tooltips
library(plotly)

ggplotly(
    umap_df %>%
        ggplot(aes(x = umap_1, y = umap_2, label = peak_id, text = peak_ann)) +
        geom_point(alpha = 0.7, colour = "black") +
        geom_text() +
        theme_classic(),
    tooltip = "text")


# MID plots
selected_exp <- c("gly", "lys", "met", "ser")
max_fraction <- 0.3

mids_tmlg <- get_mid_matrix(hmec_mi_data, "5665", selected_exp)
plot_mid_matrix(c13correct_cols(mids_tmlg)[-1, ])

write_image(
    file.path(data_path, "plots", "5665-tmlg-selected.png"),
    pmin(reverse_rows(c13correct_cols(mids_tmlg)[-1, ]), max_fraction) / max_fraction,
    colorRamp(c("white", "red"))
)

# 6269 trimethyllysine +H (known)
mids_tmlys <- get_mid_matrix(hmec_mi_data, "6269", selected_exp)
plot_mid_matrix(c13correct_cols(mids_tmlys)[-1, ])

write_image(
    file.path(data_path, "plots", "6269-tmlys-selected.png"),
    pmin(reverse_rows(c13correct_cols(mids_tmlys)[-1, ]), max_fraction) / max_fraction,
    colorRamp(c("white", "red"))
)

# 2531 glycine (known)
mids_gly <- get_mid_matrix(hmec_mi_data, "2531", selected_exp)
plot_mid_matrix(c13correct_cols(mids_gly)[-1, ])

write_image(
    file.path(data_path, "plots", "2531-gly-selected-for-tmlg.png"),
    pmin(reverse_rows(c13correct_cols(mids_gly)[-1, ]), max_fraction) / max_fraction,
    colorRamp(c("white", "red"))
)

# convolution MID tmlys * gly
mids_tmlys_gly_conv <- convolute_all(mids_tmlys, mids_gly)
plot_mid_matrix(c13correct_cols(mids_tmlys_gly_conv)[-1, ])

write_image(
    file.path(data_path, "plots", "conv_tmlys_gly.png"),
    pmin(reverse_rows(c13correct_cols(mids_tmlys_gly_conv)[-1, ]), max_fraction) / max_fraction,
    colorRamp(c("white", "red"))
)

#
# 5872 histidine compound
#

peak_skeleton %>% filter(peak_id == 5872)

find_nearest(dm, "5872", known_peak_ids, 10) %>%
    inner_join(
        peak_annotations %>% select(peak_id, known_met_id),
        by = "peak_id"
    )

hmec_peak_list[conv_index["1966", "5872"], ]

#
# 
#

