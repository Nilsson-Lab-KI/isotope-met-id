#
#  Compute the "gold standard" matrix from the fraction-derived matrix
#

source("common.R")

# fraction of carbon derived (asymmetrics)
fraction_derived <- as.matrix(
    read.table(
        file.path(input_data_path, "fraction-derived-matrix.tsv"),
        header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
    )
)
stopifnot(
    all(rownames(fraction_derived) == colnames(fraction_derived))
)

# symmetrize
fraction_derived <- pmax(fraction_derived, t(fraction_derived))
stopifnot(
    isSymmetric.matrix(fraction_derived)
)

saveRDS(
    fraction_derived,
    file.path(gold_standard_path, 'fraction_derived.rds')
)


# threshold fraction derived matrix to get a binary gold standard matrix
gold_standard <- ifelse(fraction_derived > 0.5, 1, 0)
diag(gold_standard) <- 0
stopifnot(
    all(!is.na(gold_standard))
)

saveRDS(
    gold_standard,
    file.path(gold_standard_path, 'gold_standard.rds')
)

# path lengths matrix, symmmetrized
path_lengths <- as.matrix(
    read.table(
        file.path(input_data_path, "path-length-matrix.tsv"),
        header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
    )
)

saveRDS(
    path_lengths,
    file.path(gold_standard_path, 'path_lengths.rds')
)
