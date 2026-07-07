#!/usr/bin/env Rscript
#
# Runs eem_search() on the eemCpp sample dataset (originally
# EEMcpp/accepted/eemCpp/data/{test.tab,test.gmt}) and prints the top
# expression modules, sorted by decreasing p-value (-log10 scale).
#
# Usage:
#   Rscript run_test.R
#
library(eemr)

expr_file <- system.file("extdata", "test.tab", package = "eemr")
gmt_file  <- system.file("extdata", "test.gmt", package = "eemr")

expr <- as.matrix(read.table(expr_file, header = TRUE, row.names = 1,
                              sep = "\t", check.names = FALSE))
geneSets <- read_gmt(gmt_file)

cat(sprintf("expr: %d genes x %d samples\n", nrow(expr), ncol(expr)))
cat(sprintf("geneSets: %d sets\n", length(geneSets)))

res <- eem_search(expr, geneSets)

cat("\nTop 10 expression modules:\n")
print(head(res[, c("id", "pvalue", "nSeedGenes", "nModuleGenes",
                    "relativeRadius", "absoluteRadius")], 10))

# module x sample activity matrix: mean (row-standardized) expression of each
# module's genes in each sample. By default only modules with pvalue >= 6
# (-log10 scale, i.e. p < 1e-6) are included.
activity <- module_activity(expr, res)
cat(sprintf("\nmodule activity matrix (pvalueCutoff=6): %d modules x %d samples\n",
            nrow(activity), ncol(activity)))
print(activity[, 1:5])

# no filtering: one row per module in res
activity_all <- module_activity(expr, res, pvalueCutoff = NULL)
cat(sprintf("\nunfiltered activity matrix: %d modules x %d samples\n",
            nrow(activity_all), ncol(activity_all)))
