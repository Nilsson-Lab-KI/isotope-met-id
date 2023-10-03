#
#  Code to reproduce plots in Figure 1
#

library(dplyr)
library(igraph)
library(umap)
library(ggplot2)

library(remn)

# File names
# TODO: we need to decide where to permanently store the data set

data_path <- file.path(
    "C:", "Users", "rolnil", "OneDrive - KI.SE", "Lab common", "projects",
    "reverse-engineering-deniz", "data-for-publication"
)
peak_list_path <- file.path(data_path, "nodelist-netid-pos-neg.csv")
stopifnot(file.exists(peak_list_path))
hmec_path <- file.path(data_path, "hmec-peak-areas.tsv")
stopifnot(file.exists(hmec_path))
qc_list_path <- file.path(data_path, "peak-qc-list.tsv")
stopifnot(file.exists(qc_list_path))

# List of peaks passing QC
peak_qc_list <- read.table(
    qc_list_path,
    sep = "\t",
    header = TRUE
) %>% filter(qc == 1) %>% select(peak_id)

# Read peak list
{
    hmec_peak_list <- read.delim(
        peak_list_path,
        header = T, sep = ","
    )
    # peak IDs
    hmec_peak_list <- hmec_peak_list %>% mutate(peak_id = row_number())
    # subset to QC peaks
    hmec_peak_list <- hmec_peak_list %>%
        inner_join(peak_qc_list, by = join_by(peak_id))
    # drop unused fields
    hmec_peak_list <- hmec_peak_list %>%
        select(peak_id, mass, path, ion_mode) %>%
        rename(netid_annotation = path)
    nrow(hmec_peak_list)
}


# Read peak area data from the HMEC isotope tracing experiment
hmec_peak_areas <- read.table(
    hmec_path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE
)
# Create MIData object
hmec_mi_data <- MIData(hmec_peak_areas, exp_columns = 3:length(hmec_peak_areas))
# subset the MIData object
{
    hmec_mi_data <- midata_subset(
        hmec_mi_data,
        match(
            hmec_peak_list$peak_id,
            hmec_mi_data$peak_ids,
        )
    )
}
n_peaks <- length(hmec_mi_data$peak_ids)
n_experients <- length(hmec_mi_data$experiments)

# verify that MIData peak IDs matches peak list
stopifnot(all(hmec_mi_data$peak_ids == hmec_peak_list$peak_id))


#
#  Fig 1a
#

exp_index <- match("glc", hmec_mi_data$experiments)

# peak 688 ATP-H
{
    peak_index <- get_peak_index(hmec_mi_data, "688")
    g6p_mids <- get_mids(hmec_mi_data, peak_index, exp_index)
    plot_mid_barchart(c13correct_cols(g6p_mids))
}

# peak 836 ADP-H
{
    peak_index <- get_peak_index(hmec_mi_data, "836")
    g6p_mids <- get_mids(hmec_mi_data, peak_index, exp_index)
    plot_mid_barchart(c13correct_cols(g6p_mids))
}


#
# Fig 1b
#

exp_index <- match("glc", hmec_mi_data$experiments)

# peak 1501 glucose-6P-H
{
    peak_index <- get_peak_index(hmec_mi_data, "1501")
    g6p_mids <- get_mids(hmec_mi_data, peak_index, exp_index)
    plot_mid_barchart(c13correct_cols(g6p_mids))
}

# peak 4430 UDP+H
{
    peak_index <- get_peak_index(hmec_mi_data, "4430")
    udp_mids <- get_mids(hmec_mi_data, peak_index, exp_index)
    plot_mid_barchart(c13correct_cols(udp_mids))
}

# convolution glucose-6P * UDP
{
    g6p_udp_mids <- sapply(
        1:3, function(i) convolute(g6p_mids[, i], udp_mids[, i])
    )
    plot_mid_barchart(c13correct_cols(g6p_udp_mids))
}

# peak 597 UDP-glucose-H
{
    peak_index <- get_peak_index(hmec_mi_data, "597")
    udpglc_mids <- get_mids(hmec_mi_data, peak_index, exp_index)
    plot_mid_barchart(c13correct_cols(udpglc_mids))
}

#
# Suppl Fig 1X 13C enrichment
#

glc_13c_enrichment <- sapply(
    1:n_peaks,
    function(i) isotopic_enrichment(get_avg_mid(hmec_mi_data, i, exp_index))
)
plot(sort(glc_13c_enrichment), ylim = c(0, 1))

#
# Fig 1c
#

# calculate all distances vs. UDP-glucose
# TODO: this should be a function in the remn package
{
    udpglc_index <- get_peak_index(hmec_mi_data, "597")
    dist <- rep(0.0, n_peaks)
    convolutant_index <- rep(0, n_peaks)
    for(i in 1:n_peaks) {
        assign_list[dist[i], convolutant_index[i]] <- conv_reduce(
            hmec_mi_data, i, udpglc_index, exp_index,
            f = euclidean_dist, g = which.min
        )
    }
}
# keep only distances to smaller metabolites
dist <- ifelse(
    hmec_mi_data$peak_n_atoms <= hmec_mi_data$peak_n_atoms[udpglc_index],
    dist, NA)
sum(!is.na(dist))
# distribution
hist(dist, breaks = 50)
plot(
    sort(pmin(dist, 0.5)),
    ylim = c(0.5, 0)
)

# top 20 plot
top_20_index <- order(dist)[(1:20)+1]
plot(
    dist[top_20_index],
    ylim = c(dist[last(top_20_index)], 0)
)
# corresponding peak IDs
hmec_peak_list[top_20_index, ]
hmec_mi_data$peak_n_atoms[top_20_index]
# convolutant metabolite IDs
hmec_peak_list[convolutant_index[top_20_index], ]

#
# Full distance matrix
#

# compute all pairwise distances
assign_list[dist, convolutant_index] <- conv_reduce_all(
    hmec_mi_data, exp_index,
    f = euclidean_dist, g = which.min
)

# total number of unique pairs
n_peaks*(n_peaks - 1)

# fraction of NAs (no convolutions possible)
sum(is.na(dist[lower.tri(dist)]))/ (n_peaks*(n_peaks - 1))

# distance matrix with NAs replace with largest distance
dist_no_na <- dist
dist_no_na[which(is.na(dist_no_na))] <- max_dist

# histogram of unique pairs
hist(dist[lower.tri(dist)], breaks = 100)


# cluster tree -- this looks strange
clust <- hclust(as.dist(dist_no_na), method = "average")

# write tree to pdf
#pdf(file = file.path(), height = 100, width = 200)
plot(
    as.dendrogram(clust),
    labels = 
)
#dev.off()

hmec_peak_list %>% select(peak_id, netid_annotation) %>%
    

# UMap projection, 
max_dist <- max(dist, na.rm = T)
{
    umap_proj <- umap(dist_no_na)
    umap_df <- data.frame(
        peak_id = rownames(umap_proj$layout),
        umap_1 = umap_proj$layout[,1],
        umap_2 = umap_proj$layout[,2]
    )
}
umap_df %>%
    ggplot(aes(x = umap_1, y = umap_2,
               label = peak_id,
               text = peak_id)) +
    geom_point(alpha = 0.7, colour = "black") +
    geom_text() +
    labs(x = "umap_1", y = "umap_2") +
    theme_classic()


# 
# Graph views
#
# I have not been able to make this look good.
# The "transitive closure" property of the MID distance is a problem,
# since in any pathway A + X -> B,  B + Y -> C where A, B, C are all
# similar, we get not only edges A--B and B--C but also A--C. A possible
# solution would be to discard edges like A--C if there is a path from
# A to C in the graph such where the carbon differences in each step
# is smaller than the total difference. (Here m(X) and m(Y) < m(X+Y).
# This is a "minimal equivalent graph" problem; unfortunately these are
# are generally computationally hard. A minimum spanning tree of each
# connected component might be a good alternative.
# 

# index of pairs below distance cutoff
# NOTE: the number of pairs increases sharply with distance!
dist_cutoff <- 0.1
pair_index <- which(!is.na(dist) & (dist < dist_cutoff), arr.ind = TRUE)

# make table of pairs
dist_table <- as.data.frame(pair_index) %>%
    rename(index_1 = row, index_2 = col) %>%
    mutate(
        peak_id_1 = hmec_mi_data$peak_ids[index_1],
        peak_id_2 = hmec_mi_data$peak_ids[index_2],
        distance = dist[pair_index],
        n_carbon_1 = hmec_mi_data$peak_n_atoms[index_1],
        n_carbon_2 = hmec_mi_data$peak_n_atoms[index_2],
        enrichment_1 = glc_13c_enrichment[index_1],
        enrichment_2 = glc_13c_enrichment[index_2]
        ) %>%
    filter(peak_id_1 != peak_id_2)
nrow(dist_table)  # 58,702

# direct edges so that peak 2 has more carbons
dist_table <- dist_table %>%
    filter(n_carbon_1 <= n_carbon_2)
nrow(dist_table)
# 31,147

# consider only edges where both peaks > min 13C enrichment
enrichment_cutoff = 0.2
dist_table <- dist_table %>%
    filter(pmin(enrichment_1, enrichment_2) > enrichment_cutoff)
nrow(dist_table)  # 1,823 !

# example, edge table
dist_table %>% filter(peak_id_2 == 597)

# number of incoming edges per peak
in_degree <- dist_table %>%
    select(index_2) %>%
    group_by(index_2) %>% count()
stopifnot(sum(in_degree$n) == nrow(dist_table))

# large metabolites have many incoming edges
plot(hmec_mi_data$peak_n_atoms[in_degree$index_2], in_degree$n)

# make an igraph for visualizing neighborhoods
dist_graph <- graph_from_data_frame(
    dist_table %>% select(peak_id_1, peak_id_2, distance)
)
dist_graph <- set_vertex_attr(
    dist_graph, "netid_annotation",
    value = hmec_peak_list$netid_annotation[
        match(as.integer(as_ids(V(dist_graph))), hmec_peak_list$peak_id)]
)

# create subgraph around UDP-glucose
my_subgraph <- subgraph(
    dist_graph,
    unlist(neighborhood(dist_graph, order = 2, "597", mode = "in"))
)
length(V(my_subgraph))

# igraph visualization -- not very helpful
plot.igraph(my_subgraph,
            vertex.label = hmec_peak_list$netid_annotation[hmec_peak_list$peak_id %in% subgraph_peak_ids])

# export as GraphML for visualization in yEd
write_graph(
    my_subgraph,
    file.path(data_path, "glc-597-neighborhood-2.graphml"),
    format = "graphml"
)


# export full network for cytoscape
# NOTE: this usually gives a hairball, probably more useful to look at subgraphs ...
write.table(
    hmec_peak_list,
    file.path(data_path, "glc-network-nodes.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
    dist_table,
    file.path(data_path, "glc-network-edges.tsv"),
    sep = "\t", row.names = FALSE
)

#
# Scratch
#

#  peak 2771 lipid * 556 = 597 (UDP-glucose), d(597, 2771) = 0.044
dist[
    get_peak_index(hmec_enriched_mi_data, "597"),
    get_peak_index(hmec_enriched_mi_data, "2771")
]
convolutant_index[
    get_peak_index(hmec_enriched_mi_data, "597"),
    get_peak_index(hmec_enriched_mi_data, "2771")
]
hmec_enriched_mi_data$peak_ids[294]
hmec_enriched_mi_data$peak_n_atoms[294]
# d(597, 556)
dist[
    get_peak_index(hmec_enriched_mi_data, "597"),
    get_peak_index(hmec_enriched_mi_data, "556")
]


plot_mid_barchart(get_mids(hmec_mi_data, top_20_index[2], exp_index))
plot_mid_barchart(get_mids(hmec_mi_data, convolutant_index[top_20_index[2]], exp_index))
