#
#  Code to reproduce plots in Figure 1
#

library(remn)

# Read peak area data from the HMEC isotope tracing experiment
# TODO: we need to decide where to permanently store the data set

hmec_path <- file.path(
    "C:", "Users", "rolnil", "OneDrive - KI.SE", "Lab common", "projects",
    "reverse-engineering-deniz", "data-for-publication", "hmec-peak-areas.tsv"
)
stopifnot(file.exists(hmec_path))

hmec_peak_areas <- read.table(
    hmec_path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE
)
length(hmec_peak_areas)

# Create MIData object
hmec_mi_data <- MIData(hmec_peak_areas, exp_columns = 3:length(hmec_peak_areas))
n_peaks <- length(hmec_mi_data$peak_ids)
n_experients <- length(hmec_mi_data$experiments)

# Remove false isotopomers
#hmec_mi_data_cens <- censor_false_mi(hmec_mi_data, min_experiments =  n_experients / 2)

#
#  Fig 1a
#

exp_index <- match("D-Glucose", hmec_mi_data$experiments)

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

exp_index <- match("D-Glucose", hmec_mi_data$experiments)

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
# Fig 1c
#

# calculate all distances vs. UDP-glucose
{
    udpglc_index <- get_peak_index(hmec_mi_data, "597")
    dist <- rep(0.0, n_peaks)
    middle_index <- rep(0, n_peaks)
    for(i in 1:n_peaks) {
        assign_list[dist[i], middle_index[i]] <- conv_reduce(
            hmec_mi_data, i, udpglc_index, exp_index,
            f = euclidean_dist, g = which.min
        )
    }
}
hist(dist, 50)        # full distribution

plot(sort(dist)[1:20])  # top 20

neighbors <- order(dist)[1:10]
# peak IDs for top 10
hmec_mi_data$peak_ids[neighbors]
hmec_mi_data$peak_n_atoms[neighbors]
# middle metabolite IDs for top 10
hmec_mi_data$peak_ids[middle_index[neighbors]]
    
    