#
#  Code to reproduce plots in Figure 1
#

source("common.R")


hmec_mi_data <- readRDS(file.path(mi_data_path, 'hmec_mi_data.rds'))
hmec_mi_data_censored <- readRDS(file.path(mi_data_path, 'hmec_mi_data_censored.rds'))

n_peaks <- length(hmec_mi_data$peak_ids)

plot_mid_barchart_ind_points <- function(mids)
{
   mids_df <- mids %>% t() %>% as.data.frame()
   colnames(mids_df) <- as.character(0:(ncol(mids_df)-1))

   df_long <- mids_df %>%
      mutate(replicate = 1:3) %>%
      pivot_longer(cols = -replicate,
                   names_to = "mi",
                   values_to = "fraction") %>%
      mutate(mi = factor(mi, levels = unique(mi)))  # avoid re-ordering labels

   # Calculate means for bars
   means <- df_long %>%
      group_by(mi) %>%
      summarise(mean_fraction = mean(fraction))
   
   # Plot
   ggplot() +
      geom_col(data = means, aes(x = mi, y = mean_fraction)) +
      geom_point(
         data = df_long, aes(x = mi, y = fraction), 
         position = position_dodge2(width = 0.7),
         size = 10, color = "black") +
      theme_classic() +
      labs(y = "MI fraction")
}


#
#  Figure 1a
#

# peak 1917 citrate-H
{
    mids <- get_mids(hmec_mi_data, "1917", "glc")
    plot_mid_barchart_ind_points(c13correct_cols(mids))
}

# peak 2064 aconitate-H (candidate)
{
   mids <- get_mids(hmec_mi_data, "2064", "glc")
   plot_mid_barchart_ind_points(c13correct_cols(mids))
}


#
# Figure 1b
#

# peak 1501 glucose-6P-H
{
   g6p_mids <- get_mids(hmec_mi_data, "1501", "glc")
   plot_mid_barchart_ind_points(c13correct_cols(g6p_mids))
}

# peak 4430 UDP+H
{
   udp_mids <- get_mids(hmec_mi_data, "4430", "glc")
   plot_mid_barchart_ind_points(c13correct_cols(udp_mids))
}


# convolution glucose-6P * UDP
{
   g6p_udp_mids <- convolute_all(g6p_mids, udp_mids)
   plot_mid_barchart_ind_points(c13correct_cols(g6p_udp_mids))
}


# peak 597 UDP-glucose-H
{
   udpglc_mids <- get_mids(hmec_mi_data, "597", "glc")
   plot_mid_barchart_ind_points(c13correct_cols(udpglc_mids))
}



#
# ED Figure 1a,b
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
# Figure 1c, toy example
#

mid_A <- c(0.7, 0.02, 0.08, 0.2)
mid_B <- c(0.35, 0.05, 0.4, 0.1, 0.01, 0.09)
mid_C1 <- c(0.9, 0.1, 0)
mid_C2 <- c(0.5, 0.05, 0.45)
mid_C3 <- c(0.1, 0.05, 0.85)

barplot(mid_A, ylim = c(0,1))

barplot(mid_C1, ylim = c(0,1))
barplot(mid_C2, ylim = c(0,1))
barplot(mid_C3, ylim = c(0,1))

barplot(
   matrix(c(convolute(mid_A, mid_C1), mid_B), nrow = 2, byrow = TRUE),
   ylim = c(0,1), beside = TRUE)
barplot(
   matrix(c(convolute(mid_A, mid_C2), mid_B), nrow = 2, byrow = TRUE),
   ylim = c(0,1), beside = TRUE)
barplot(
   matrix(c(convolute(mid_A, mid_C3), mid_B), nrow = 2, byrow = TRUE),
   ylim = c(0,1), beside = TRUE)

euclidean_dist(convolute(mid_A, mid_C1), mid_B)
euclidean_dist(convolute(mid_A, mid_C2), mid_B)
euclidean_dist(convolute(mid_A, mid_C3), mid_B)


#
# Figure 1d
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

# histogram, clip at >0.5
hist(pmin(udpglc_dist, 0.5), breaks = 50)

# plot, clip at >0.5
plot(
    sort(pmin(udpglc_dist, 0.5)),
    ylim = c(0.5, 0)
)

#
# Figure 1e
#

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

# top peaks and convolutants
udpglc_top_neighbors <- data.frame(
    peak_id = hmec_mi_data$peak_ids[udpglc_top_index]
)

udpglc_top_convolutants <- data.frame(
    peak_id = hmec_mi_data$peak_ids[convolutant_index[udpglc_top_index]]
)

# candidate annotations from HMDB
peak_hmdb_compound <- read_peak_hmdb_compound()

udpglc_top_neighbors %>%
    inner_join(
        peak_hmdb_compound,
        by = "peak_id"
    )

udpglc_top_convolutants %>%
    inner_join(
        peak_hmdb_compound,
        by = "peak_id",
        relationship = "many-to-many"
    )

