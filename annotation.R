
library(dplyr)
library(rcdk)
library(stringr)

library(remn)

source("common.R")

#
# Setup
#

{
    project_path <-  file.path(
        "C:", "Users", "rolnil", "OneDrive - KI.SE", "Lab common", "projects",
        "reverse-engineering-deniz"
    )
    data_path <- file.path(project_path, "data-for-publication")
    peak_list_path <- file.path(data_path, "nodelist-netid-pos-neg.csv")
    stopifnot(file.exists(peak_list_path))
    hmec_path <- file.path(data_path, "hmec-peak-areas.tsv")
    stopifnot(file.exists(hmec_path))
    qc_list_path <- file.path(data_path, "peak-qc-list.tsv")
    stopifnot(file.exists(qc_list_path))
}

{
    proton_mass <- 1.007276466621
    h2o_mass <- 18.010564683
    nh4_mass <- 18.034374132
}

#
# Structure database
#

# Read list of compounds exported from HMDB
# This contains distinct inchi keys, and neutral masses in all cases
structure <- read.table(
    file.path(project_path, "hmdb", "hmdb_compounds.tsv"),
    header = TRUE, sep = "\t", quote = "", comment.char = ""
)
structure %>% distinct(inchi_key) %>% nrow() == structure %>% nrow()

# filter by mass to the measured mass range
structure <- structure %>% filter(neutral_mass > 50 & neutral_mass < 1000)
# drop structures without SMILES
structure <- structure %>% filter(smiles != "")
structure%>% count

# unique molecular skeletons
structure <- structure %>%
    mutate(skeleton = str_sub(inchi_key, end = 14))

structure %>% distinct(skeleton) %>% nrow
structure %>% distinct(smiles) %>% nrow

# verify that all skeletons have the same neutral mass, up to rounding errors
structure %>% 
    mutate(neutral_mass = round(neutral_mass, digits = 5)) %>%
    distinct(skeleton, neutral_mass) %>%
    group_by(skeleton) %>% filter(n() > 1) %>% nrow == 0

# table of unique skeletons, with a representative SMILES and mass
skeleton <- structure %>%
    select(skeleton, smiles, neutral_mass) %>%
    distinct(skeleton, .keep_all = TRUE)

# Alternative: structure database from KEGG / KMRN
# NOTE: not all compounds here are neutral, and in such cases
# the masses given are not the neutral mass!
# structure <- read.table(
#     "X:/datasets/kgmn-zhou-2022/KMRN_known/node_table.tsv",
#     header = TRUE, sep = "\t", quote = "", comment.char = ""
# )

# drop some duplicates
# structure <- structure %>% distinct(inchi_key, .keep_all = TRUE)

# check charge field, drop non-neutral structures (everything except "N")
# structure %>% select(inchikey) %>%
#     mutate(charge_key = str_sub(structure$inchikey, -1)) %>%
#     group_by(charge_key) %>% count()
# 
# structure <- structure %>% filter(str_sub(structure$inchikey, -1) == "N")

# number of unique structures
# structure %>% select(inchikey1, cpd_name) %>% distinct() %>% count()
# number of unique molecular skeletons
# structure %>% select(inchikey1) %>% distinct() %>% count()

# collapse to unique molecular skeleton (ignoring charge and stereochemistry)
# skeleton <- structure %>%
#     select(inchikey, inchikey1, smiles, cpd_name, monoisotopic_mass) %>% distinct() %>%
#     group_by(inchikey1) %>% summarize(
#         candidate_names = str_c(cpd_name, collapse = "|"),
#         monoisotopic_mass = first(monoisotopic_mass),
#         repr_key = first(inchikey),
#         repr_smiles = first(smiles)
#     )


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
    # recover the m/z from neutral mass
    hmec_peak_list <- hmec_peak_list %>%
        mutate(mz = mass + proton_mass*ion_mode)
    # drop unused fields
    hmec_peak_list <- hmec_peak_list %>%
        select(peak_id, mass, path, ion_mode, mz) %>%
        rename(netid_annotation = path)
    nrow(hmec_peak_list)
}

# map each peak to possible ions formed
# NOTE: not all possible in each ion mode
ion_modification <- data.frame(
    mod_id = c("+H", "-H", "+NH4", "+H-H2O", "-H-H2O"),
    mod_mass = c(proton_mass, -proton_mass, nh4_mass, proton_mass - h2o_mass, -proton_mass - h2o_mass),
    ion_mode = c(1, -1, 1, 1, -1)
)

ppm_tolerance <- 10
hmec_peak_ions <- hmec_peak_list %>%
    inner_join(
        ion_modification,
        by = join_by(ion_mode),
        relationship = "many-to-many"
    ) %>%
    mutate(
        mass_min = (mz - mod_mass) * (1 - ppm_tolerance*1e-6),
        mass_max = (mz - mod_mass) * (1 + ppm_tolerance*1e-6)
    )

# map ions to candidate molecular skeletons
peak_skeleton <- skeleton %>%
    inner_join(
        hmec_peak_ions,
        by = join_by(between(neutral_mass, mass_min, mass_max))
    )

# number of peaks with a matching candidate structure
# NOTE: still have a number of in source fragments in here

# number of candidate skeletons per peak
peak_skeleton %>%
    group_by(peak_id) %>% count(name = "n_candidates") %>%
    {.$n_candidates} %>% sort %>% plot

peak_skeleton %>%
    group_by(peak_id) %>% count(name = "n_candidates") %>%
    group_by(n_candidates) %>% count %>% plot

# total number of skeletons to consider
peak_skeleton %>% distinct(skeleton) %>% nrow

# list of known compounds, with inchi keys
known_compounds <- read.table(
    file.path(data_path, "qc721-known-compounds.tsv"),
    header = TRUE, sep = "\t", quote = "", na.strings = ""
) %>%
    filter(!is.na(met_id)) %>%
    mutate(skeleton = str_sub(inchi_key, 1, 14))

known_compounds %>% nrow

# check for known inchi keys missing in the database
known_compounds %>% anti_join(structure, by = "inchi_key") %>% nrow
# known_compounds <- known_compounds %>% filter(inchi_key %in% structure$inchikey)

# number of unique known skeletons
known_compounds %>% distinct(skeleton) %>% nrow

# create CDK structure objects for each candidate skeleton
cdk_mols <- parse.smiles(
    peak_skeleton %>% distinct(skeleton, .keep_all = TRUE) %>% {.$smiles}
)
names(cdk_mols) <- peak_skeleton %>% distinct(skeleton) %>% {.$skeleton}


#
# read HMEC data and make distance matrix
#

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

stopifnot(all(hmec_mi_data$peak_ids == hmec_peak_list$peak_id))

# distance matrix based on all experiments except unlabeled
# this takes a few minutes
# assign_list[dm, conv_index] <- conv_reduce_all(
#     hmec_mi_data, 1:(n_experients - 1),
#     f = remn::euclidean_sum_dist,
#     g = which.min,
#     impute = TRUE
# )

# save this object to speed up future work
# save(dm, conv_index, file = file.path(data_path, "dm-conv_index.RData"))

# read precomputed matrix
load(file.path(data_path, "dm-conv_index.RData"))

#
# annotate peaks
#

# 1. for each unknown peak, find all known peaks within a distance radius
# and their convolutants
find_known <- function(dm, conv_index, known_ids, max_distance = 1.0)
{
    unknown_ids <- setdiff(colnames(dm), as.character(known_ids))
    known_dm <- dm[as.character(known_ids), unknown_ids]
    known_conv_index <- conv_index[as.character(known_ids), unknown_ids]
    pair_index <- which(known_dm < max_distance, arr.ind = TRUE)
    data.frame(
        known_id = as.integer(rownames(known_dm)[pair_index[, 1]]),
        unknown_id = as.integer(colnames(known_dm)[pair_index[, 2]]),
        conv_id = as.integer(colnames(known_dm)[conv_index[pair_index]])
    )
}

# table with all matching knowns peaks for each unknown, and associated skeletons
peak_candidates <- find_known(dm, conv_index, known_compounds$peak_id) %>%
    inner_join(
        known_compounds %>% select(peak_id, skeleton, met_id),
        by = join_by(known_id == peak_id)
    ) %>%
    rename(known_skeleton = skeleton, known_met_id = met_id) %>%
    inner_join(
        peak_skeleton %>% select(peak_id, skeleton),
        by = join_by(unknown_id == peak_id),
        relationship = "many-to-many"
    ) %>%
    rename(candidate_skeleton = skeleton) %>%
    select(
        unknown_id, candidate_skeleton, known_id, known_skeleton, known_met_id)

# number of unknown peaks associated with at least one known
peak_candidates %>% distinct(unknown_id) %>% nrow

# distribution of no. of similar known peaks per unknown peak
peak_candidates %>% select(unknown_id, known_id) %>% distinct() %>%
    group_by(unknown_id) %>% count(name = 'n_known') %>%
    group_by(n_known) %>% count() %>% plot

# total number of candidate structures
peak_candidates %>% distinct(candidate_skeleton) %>% nrow()

# number of structure pairs to consider for MCS
peak_candidates %>% distinct(known_skeleton, candidate_skeleton) %>% nrow()


# 2. for each candidate structure we add the MCS similarity

# number of atoms in MCS normalized to maximum
mcs_similarity <- Vectorize(
    function(skeleton_1, skeleton_2)
    {
        mol_1 <- cdk_mols[[skeleton_1]]
        mol_2 <- cdk_mols[[skeleton_2]]
        return(
            get.atom.count(get.mcs(mol_1, mol_2)) /
                max(get.atom.count(mol_1), get.atom.count(mol_2))
        )
    }
)

# this step is slow due to MCS computation
# NOTE: this recomputes MCS when multiple peaks map to the same structures,
# could be optimized
peak_candidates <- peak_candidates %>%
    mutate(
        mcs_similarity = mcs_similarity(candidate_skeleton, known_skeleton)
    )

# keep structures above some minimal MCS similarity and
# for each unknown, sekect candidate structures with best MCS
min_mcs = 0.4
selected_candidates <- peak_candidates %>%
    filter(mcs_similarity >= min_mcs) %>%
    inner_join(
        peak_candidates %>%
            group_by(unknown_id) %>% summarise(max_sim = max(mcs_similarity)),
        by = join_by(unknown_id == unknown_id, mcs_similarity == max_sim)
    )


# number of unknowns with at least one prediction
selected_candidates %>% distinct(unknown_id) %>% count

# annotate candidate list
annotated_candidates <- selected_candidates %>%
    select(unknown_id, candidate_skeleton) %>%
    inner_join(
        structure %>% select(skeleton, name),
        by = join_by(candidate_skeleton == skeleton),
        relationship = "many-to-many"
    ) %>%
    group_by(unknown_id, candidate_skeleton) %>%
    summarize(names = str_c(name, collapse = "|"))


# look up specific annotations

peak_candidates %>% filter(unknown_id == 2031) %>%
    inner_join(
        structure %>% select(skeleton, name),
        by = join_by(candidate_skeleton == skeleton),
        relationship = "many-to-many"
    )


selected_candidates %>% filter(unknown_id == 2031) %>%
    inner_join(
        structure %>% select(skeleton, name),
        by = join_by(candidate_skeleton == skeleton),
        relationship = "many-to-many"
    )


