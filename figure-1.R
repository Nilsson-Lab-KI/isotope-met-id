#
#  Code to reproduce plots in Figure 1
#

library(dplyr)

library(remn)

# File names
# TODO: we need to decide where to permanently store the data set

data_path <- file.path(
    "C:", "Users", "rolnil", "OneDrive - KI.SE", "Lab common", "projects",
    "reverse-engineering-deniz", "data-for-publication"
)
hmec_path <- file.path(data_path, "hmec-peak-areas.tsv")
stopifnot(file.exists(hmec_path))
qc_list_path <- file.path(data_path, "peak-qc-list.tsv")
stopifnot(file.exists(qc_list_path))

# Read peak area data from the HMEC isotope tracing experiment
hmec_peak_areas <- read.table(
    hmec_path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE
)

# Create MIData object
hmec_mi_data <- MIData(hmec_peak_areas, exp_columns = 3:length(hmec_peak_areas))
n_peaks <- length(hmec_mi_data$peak_ids)
n_experients <- length(hmec_mi_data$experiments)

# TODO: Remove false isotopomers
# this should probably be a separate a processing script 
# manual QC will of course have to be just a list of peak IDs
#hmec_mi_data_cens <- censor_false_mi(hmec_mi_data, min_experiments =  n_experients / 2)

# Subset to final QC'ed peak list
peak_qc_list <- read.table(
    qc_list_path,
    sep = "\t",
    header = TRUE
)

# subset the MIData object
{
    qc_peak_ids <- peak_qc_list %>% filter(qc == 1) %>% select(peak_id)
    hmec_qc_mi_data <- midata_subset(
        hmec_mi_data,
        match(
            qc_peak_ids[, 1],
            hmec_mi_data$peak_ids,
        )
    )
    n_qc_peaks <- length(hmec_qc_mi_data$peak_ids)
}

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
# TODO: this could be a function in the package
{
    udpglc_index <- get_peak_index(hmec_qc_mi_data, "597")
    dist <- rep(0.0, n_qc_peaks)
    convolutant_index <- rep(0, n_qc_peaks)
    for(i in 1:n_qc_peaks) {
        assign_list[dist[i], convolutant_index[i]] <- conv_reduce(
            hmec_qc_mi_data, i, udpglc_index, exp_index,
            f = euclidean_dist, g = which.min
        )
    }
}
# full distribution
hist(dist, breaks = 30)
plot(sort(dist), ylim = c(max(dist, na.rm = TRUE), 0))

# top 20 plot
top_20_index <- order(dist)[(1:20)+1]
plot(dist[top_20_index], ylim = c(dist[last(top_20_index)], 0))
# corresponding peak IDs
hmec_qc_mi_data$peak_ids[top_20_index]
hmec_qc_mi_data$peak_n_atoms[top_20_index]
# middle metabolite IDs for top 10
hmec_qc_mi_data$peak_ids[convolutant_index[top_20_index]]
    
    