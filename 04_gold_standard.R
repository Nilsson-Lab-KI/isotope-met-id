#
#  Compute the "gold standard" matrix from the fraction-derived matrix
#

source("common.R")


frac_derived <- as.matrix(
    read.table(
        file.path(input_data_path, "fraction-derived-matrix.tsv"),
        header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
    )
)
stopifnot(
    all(rownames(frac_derived) == colnames(frac_derived))
)

# symmetrize
frac_derived_symm <- pmax(frac_derived, t(frac_derived))
stopifnot(
    isSymmetric.matrix(frac_derived_symm)
)

# threshold fraction derived matrix to get a binary gold standard matrix
gold_standard <- ifelse(frac_derived > 0.5, 1, 0)
stopifnot(
    all(!is.na(gold_standard))
)

saveRDS(
    gold_standard,
    file.path(gold_standard_path, 'gold_standard.rds')
)

