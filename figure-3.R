#
#  Figure 3

source("common.R")

library(reshape2)
library(stringr)


#
# ED Figure 3b
#

fraction_derived <- readRDS(file.path(gold_standard_path, 'fraction_derived.rds'))
metabolite_pathway <- read_metabolite_pathway()

frac_derived_long <- melt(
    fraction_derived,
    varnames = c("metabolite_1", "metabolite_2"),
    value.name = "fraction_derived"
)

# pairs of metabolites in the same pathway, lower-triangular
pathway_pairs <- metabolite_pathway %>%
    rename(metabolite_1 = id) %>%
    inner_join(
        metabolite_pathway %>%
            rename(metabolite_2 = id),
        by = 'pathway',
        relationship = 'many-to-many'
    ) %>%
    filter(metabolite_1 < metabolite_2)


# fraction derived between metabolite in same pathway, ignoring compartments
frac_derived_same_pathway <- frac_derived_long %>%
    mutate(
        metabolite_1 = str_sub(metabolite_1, end = -3),
        metabolite_2 = str_sub(metabolite_2, end = -3)
    ) %>%
    inner_join(
        pathway_pairs,
        by = c("metabolite_1", "metabolite_2"),
        relationship = 'many-to-many'
    )

frac_derived_diff_pathway <- frac_derived_long %>%
    anti_join(
        frac_derived_same_pathway,
        by = c("metabolite_1", "metabolite_2")
    )
    
# plot histograms
hist(
    frac_derived_same_pathway$fraction_derived,
    breaks = 50, freq = FALSE)

hist(
    frac_derived_diff_pathway$fraction_derived,
    breaks = 50, freq = FALSE)


#
# ED Figure 3c Comparison with path lengths
#
# for computation of path lengths from atom map, see Mathematica code
#

path_lengths = readRDS(file.path(gold_standard_path, 'path_lengths.rds'))

local({
    original_par = par(no.readonly = TRUE)
    par(pty = "s")
    plot(
        x = fraction_derived,
        y = path_lengths,
        col = alpha("black", 0.05)
    )
    par(original_par)
})

