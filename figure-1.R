#
#  Code to reproduce plots in Figure 1
#


library(remn)

source("common.R")

hmec_mi_data <- readRDS(file.path(mi_data_path, 'hmec_mi_data.rds'))
n_peaks <- length(hmec_mi_data$peak_ids)

#
#  Fig 1a
#

# peak 1917 citrate-H
{
    g6p_mids <- get_mids(hmec_mi_data, "1917", "glc")
    plot_mid_barchart(c13correct_cols(g6p_mids))
}

# peak 2064 aconitate-H (candidate)
{
    g6p_mids <- get_mids(hmec_mi_data, "2064", "glc")
    plot_mid_barchart(c13correct_cols(g6p_mids))
}


#
# Fig 1b
#

# peak 1501 glucose-6P-H
{
    g6p_mids <- get_mids(hmec_mi_data, "1501", "glc")
    plot_mid_barchart(c13correct_cols(g6p_mids))
}

# peak 4430 UDP+H
{
    udp_mids <- get_mids(hmec_mi_data, "4430", "glc")
    plot_mid_barchart(c13correct_cols(udp_mids))
}


# convolution glucose-6P * UDP
{
    g6p_udp_mids <- convolute_all(g6p_mids, udp_mids)
    plot_mid_barchart(c13correct_cols(g6p_udp_mids))
}


# peak 597 UDP-glucose-H
{
    udpglc_mids <- get_mids(hmec_mi_data, "597", "glc")
    plot_mid_barchart(c13correct_cols(udpglc_mids))
}

#
# Suppl Fig 1a,b 13C enrichment
#

enrichment_vector <- function(mi_data, exp_name)
{
    sapply(
        1:length(mi_data$peak_ids),
        function(i) isotopic_enrichment(get_avg_mid(mi_data, i, exp_name))
    )
}

enrichment <- data.frame(
    glc = enrichment_vector(hmec_mi_data, "glc"),
    unlabeled = enrichment_vector(hmec_mi_data, "unlabeled"),
    row.names = hmec_mi_data$peak_ids
) %>% arrange(glc)

{
    plot(x = (1:n_peaks), y = enrichment$glc, ylim = c(0,1))
    points(x = (1:n_peaks), y = enrichment$unlabeled, col = "blue")
}

# Censor false MIs
hmec_mi_data_censored <- censor_false_mi(
    hmec_mi_data,
    min_experiments = 10
) 

enrichment_censored <- data.frame(
    glc = enrichment_vector(hmec_mi_data_censored, "glc"),
    unlabeled = enrichment_vector(hmec_mi_data_censored, "unlabeled"),
    row.names = hmec_mi_data_censored$peak_ids
) %>% arrange(glc)

{
    plot(x = (1:n_peaks), y = enrichment_censored$glc, ylim = c(0,1))
    points(x = (1:n_peaks), y = enrichment_censored$unlabeled, col = "blue")
}



#
# Fig 1c
#


# calculate all distances vs. UDP-glucose
# TODO: this should be a function in the remn package
{
    udpglc_index <- get_peak_index(hmec_mi_data_censored, "597")
    udpglc_dist <- rep(0.0, n_peaks)
    convolutant_index <- rep(0, n_peaks)
    for(i in 1:n_peaks) {
        assign_list[udpglc_dist[i], convolutant_index[i]] <- conv_reduce(
            hmec_mi_data, i, udpglc_index, "glc",
            f = euclidean_dist, g = which.min
        )
    }
}

# keep only distances to smaller metabolites
# dist <- ifelse(
#     hmec_mi_data$peak_n_atoms <= hmec_mi_data$peak_n_atoms[udpglc_index],
#     dist, NA)
# sum(!is.na(dist))
# distribution

hist(udpglc_dist, breaks = 50)

plot(
    sort(pmin(udpglc_dist, 0.5)),
    ylim = c(0.5, 0)
)

# zoom in on d < 0.1
udpglc_top_index <- intersect(
    order(udpglc_dist),
    which(udpglc_dist < 0.1)
)[-1]
length(udpglc_top_index)

plot(
    udpglc_dist[udpglc_top_index],
    ylim = c(0.1, 0)
)

# corresponding peaks
udpglc_top_neighbors <- hmec_peak_list[udpglc_top_index, ]
udpglc_top_neighbors

# convolutant peaks
top_convolutants <- hmec_peak_list[convolutant_index[udpglc_top_index], ]


