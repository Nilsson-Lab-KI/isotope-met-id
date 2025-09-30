#
# Summarize manual annotation
# NOTE: remove this files when we have moved all relevant stuff to figure-4.R
#

library(dplyr)
library(stringr)

source("common.R")

# TODO: remove this to avoid depending on readxl
#       replace with read_hmec_peak_list()
# peak_annotations <- readxl::read_excel(
#     file.path(project_path, "reverse-engineering", "Results", "hmec", "node_list_qc_tree_ordered_RN.xlsx")
# ) %>% mutate(skeleton = str_sub(inchi_key, end = 14)) 



# number of peaks annotated as lipids
n_lipids <- peak_annotations %>% filter(is_lipid == 1) %>% nrow

# no lipids are annotated with an inchi
peak_annotations %>% filter((is_lipid == 1) & !is.na(inchi_key)) %>% nrow

# total no. peaks with an MS2 spectrum
peak_annotations %>% filter(has_ms2 == 1) %>% nrow

# total no. peaks with a GNPS library hit
peak_annotations %>% filter(!is.na(gnps_library_name)) %>% nrow

# no. lipids with a GNPS library hit
peak_annotations %>%
    filter(is_lipid == 1 & !is.na(gnps_library_name)) %>% nrow

# check that we don't have NAs messing up counts
stopifnot(peak_annotations %>% filter(is.na(is_known)) %>% nrow == 0)


#
# Excluding lipids
#

# remove lipids
peak_annotations <- peak_annotations %>% filter(is_lipid == 0)

# total number of peaks annotated with an inchi
peak_annotations %>% filter(!is.na(inchi_key)) %>% nrow
# number of unique inchi
peak_annotations %>% filter(!is.na(inchi_key)) %>% 
    distinct(inchi_key) %>% nrow

# previously known peaks
peak_annotations %>% filter(is_known == 1) %>% nrow

peak_annotations %>% filter(is_known == 1) %>%
    distinct(inchi_key) %>% nrow

# new annotated peaks
peak_annotations %>% filter(is_known == 0) %>%
    filter(!is.na(inchi_key)) %>% nrow

peak_annotations %>% filter(is_known == 0) %>%
    filter(!is.na(inchi_key)) %>%
    distinct(inchi_key) %>% nrow

# number of in-source fragments
peak_annotations %>% filter(is_adduct_fragment_manual == 1) %>% nrow

# number of dipeptides
peak_annotations %>% filter(is_dipeptide == 1) %>% nrow

#
# Comparison with MS2 matching
#


# total number of peaks with an MS2 spectrum
peak_annotations %>% filter(has_ms2 == 1) %>% nrow

# total number of peaks with a GNPS library hit (lipids excluded)
peak_annotations %>% filter(!is.na(gnps_library_name)) %>% nrow

# known peaks with an MS2 spectrum
peak_annotations %>% filter(is_known == 1 & has_ms2 == 1) %>% nrow
# known peaks with a GNPS library hit
peak_annotations %>% filter(is_known == 1 & !is.na(gnps_library_name)) %>% nrow

# no. known metabolites with an MS2 spectrum
peak_annotations %>% filter(is_known == 1 & has_ms2 == 1) %>% 
    distinct(inchi_key) %>% nrow
# no. known metabolites with a GNPS hit (for at least one peak)
peak_annotations %>%
    filter(is_known == 1 & !is.na(gnps_library_name)) %>%
    distinct(inchi_key) %>% nrow

# known metabolites without a GNPS hit
peak_annotations %>%
    filter(is_known == 1) %>%
    select(inchi_key, gnps_library_name) %>%
    group_by(inchi_key) %>%
    filter(!is.na(gnps_library_name))

# no. unknown preaks
n_unknown <- peak_annotations %>% filter(is_known == 0) %>% nrow

# no. unknown peaks with an MID annotation
n_unknown_mid_ann <- peak_annotations %>%
    filter(is_known == 0 & !is.na(inchi_key)) %>% nrow

# unknown metabolites
n_unknown_ms2_ann <-peak_annotations %>%
    filter(is_known == 0 & !is.na(gnps_library_name)) %>% nrow

# bar chart of fractions
barplot(
    c(1, n_unknown_mid_ann / n_unknown, n_unknown_ms2_ann / n_unknown)
)


# no. unknown peaks with both MID annotation and GNPS library hit
peak_annotations %>%
    filter(is_known == 0 & !is.na(inchi_key) & !is.na(gnps_library_name)) %>% nrow

279/538

#
# Scratch
#

structure <- read.table(
    file.path(project_path, "hmdb", "hmdb_compounds.tsv"),
    header = TRUE, sep = "\t", quote = "", comment.char = ""
)
structure %>% distinct(inchi_key) %>% nrow() == structure %>% nrow()

# filter by mass to the measured mass range
structure <- structure %>% filter(neutral_mass > 50 & neutral_mass < 1000)
# drop structures without SMILES
structure <- structure %>% filter(smiles != "")
# add inchi key skeleton
structure <- structure%>% mutate(skeleton = str_sub(inchi_key, end = 14))

structure%>% count




# Venn diagram for non-lipids

# 538 non-lipids peaks total = 50^2
# 101 known peaks 
sqrt(50^2 * 101/538)

# 437 unknown non-lipids peaks total = 50^2
# 270 MID annotations
sqrt(50^2 * 270/437)
# 16 MS2 library hits
sqrt(50^2 * 16/437)


# MID-annotated peaks with an MS2 spectrum
peak_annotations %>%
    filter(is_known == 0 & !is.na(inchi_key) & has_ms2 == 1) %>% nrow

# MID-annotated peaks with a GNPS library hit
peak_annotations %>% filter(is_known == 0) %>%
    select(peak_id, hypothesis, is_adduct_fragment_manual, inchi_key, gnps_library_name, gnps_library_inchi) %>%
    filter(!is.na(inchi_key) & !is.na(gnps_library_name))

barplot(c(1.0, 54/183, 270/437, 16/437))
