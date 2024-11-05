#
# Preprocess input data
#

source("common.R")

library(stringr)

#
# HMEC experiment peak list
#

# peak list processed by NetID
netid_peak_list <- read.delim(
    file.path(input_data_path, "nodelist-netid-pos-neg.csv"),
    header = T, sep = ","
)
# number of putative metabolites (excluding fragments, artefacts)
netid_peak_list %>%
    filter(class %in% c('Metabolite', 'Putative metabolite')) %>% nrow()

netid_peak_list <- netid_peak_list %>%
    # define peak IDs as row number
    mutate(peak_id = as.character(row_number())) %>%
    # recover the m/z from neutral mass +/- proton
    mutate(mz = mass + proton_mass*ion_mode) %>%
    # drop unused fields
    select(peak_id, mass, ion_mode, mz)

# list of 721 peaks passing QC
qc_peak_list <- read_tsv(file.path(input_data_path, "qc-peaks-with-flags.tsv")) %>%
    mutate(peak_id = as.character(peak_id)) %>%
    mutate(known_met_id = na_if(known_met_id, "")) %>%
    mutate(inchi_key = na_if(inchi_key, "")) %>%
    mutate(skeleton = na_if(skeleton, ""))

# merge peak lists
hmec_peak_list <-  netid_peak_list %>%
    inner_join(qc_peak_list, by = join_by(peak_id))

write.table(
    hmec_peak_list, file.path(preprocessed_data_path, 'hmec_peak_list.tsv'),
    sep = '\t', row.names = FALSE, quote = FALSE
)


# list of compounds exported from HMDB
hmdb_compounds <- read_tsv(file.path(input_data_path, "hmdb_compounds.tsv"))
# verify inchi_key is unique
stopifnot(
    hmdb_compounds %>% distinct(inchi_key) %>% nrow() == hmdb_compounds %>% nrow())

# filter by mass to the measured mass range
hmdb_compounds <- hmdb_compounds %>% filter(neutral_mass > 50 & neutral_mass < 1000)
hmdb_compounds %>% count

# molecular skeletons are defined by the first 14 characters of the inchi key,
# discarding stereochemistry and charge information
hmdb_compounds <- hmdb_compounds %>%
    mutate(skeleton = str_sub(inchi_key, end = 14))

# verify that all molecular skeletons have the same neutral mass, up to rounding errors
stopifnot(
    hmdb_compounds %>% 
        mutate(neutral_mass = round(neutral_mass, digits = 5)) %>%
        distinct(skeleton, neutral_mass) %>%
        group_by(skeleton) %>% filter(n() > 1) %>% nrow == 0)

# collapse to table of unique molecular skeletons, with a representative inchi and mass
# NOTE: some compounds have multiple entries in HMDB, with different inchis
#       (e.g. with and w/o stereochemistry information)
#       this will retain only one inchi key (arbitrarily)
hmdb_compounds <- hmdb_compounds %>%
    select(skeleton, inchi_key, name, hmdb_accession, neutral_mass) %>%
    distinct(skeleton, .keep_all = TRUE)
hmdb_compounds %>% count


# map each peak to possible ions formed

# ion modifications we consider
ion_modification <- data.frame(
    mod_id = c("+H", "-H", "+NH4", "+H-H2O", "-H-H2O"),
    mod_mass = c(proton_mass, -proton_mass, nh4_mass, proton_mass - h2o_mass, -proton_mass - h2o_mass),
    ion_mode = c(1, -1, 1, 1, -1)
)

# generate possible ions from peak list
hmec_peak_ions <- hmec_peak_list %>%
    select(peak_id, mz, ion_mode) %>%
    inner_join(
        ion_modification,
        by = join_by(ion_mode),
        relationship = "many-to-many"
    ) %>%
    mutate(
        mass_min = (mz - mod_mass) * (1 - ppm_tolerance*1e-6),
        mass_max = (mz - mod_mass) * (1 + ppm_tolerance*1e-6)
    )

# map ions to candidate HMDB annotations
peak_hmdb_compound <- hmdb_compounds %>%
    inner_join(
        hmec_peak_ions,
        by = join_by(between(neutral_mass, mass_min, mass_max))) %>%
    select(
        peak_id, ion_mode, mz,
        hmdb_accession, name, inchi_key, skeleton, neutral_mass,
        mod_id, mod_mass) %>%
    arrange(peak_id)

# save peak_compound table as input data
write.table(
    peak_hmdb_compound,
    file.path(preprocessed_data_path, "peak_hmdb_compound.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

# read peak areas for all peaks
{
    hmec_peak_areas <- read_tsv(file.path(input_data_path, "hmec-peak-areas.tsv"))
    # order tracer data columns alphabetic by tracer, unlabeled data last
    tracer_data_col_names = colnames(hmec_peak_areas)[3:65]
    column_order <- c(1:2, 2 + order(tracer_data_col_names), 66:68)
    # subset to QC peaks
    hmec_peak_areas <- hmec_peak_areas[
        which(hmec_peak_areas$peak_id %in% hmec_peak_list$peak_id),
        column_order
    ]
}

# verify peak IDs match
stopifnot(
    all(
        hmec_peak_areas %>% select(peak_id) %>% distinct() ==
        qc_peak_list %>% select(peak_id)))

# write table of peak areas for QC peaks
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


