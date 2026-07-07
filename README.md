# eemr

R bindings (via Rcpp) for the original **eemCpp** coherence-based expression
module (EEM) search tool. Given an expression matrix and a collection of
candidate gene sets, `eem_search()` finds the coherently expressed subset of
each gene set ("expression module") and assigns it a significance p-value,
without leaving R (no temp files, no CLI binary).

## Install

```r
# from a local clone
devtools::install("path/to/eemr")
```

## Usage

```r
library(eemr)

expr <- as.matrix(read.table("expr.tab", header = TRUE, row.names = 1,
                              sep = "\t", check.names = FALSE))
geneSets <- read_gmt("genesets.gmt")

res <- eem_search(expr, geneSets)
head(res[, c("id", "pvalue", "nSeedGenes", "nModuleGenes")])

# module x sample activity matrix (mean row-standardized expression of each
# module's genes), restricted by default to modules with pvalue >= 6
# (-log10 scale, i.e. p < 1e-6)
activity <- module_activity(expr, res)
```

See [`inst/examples/run_test.R`](inst/examples/run_test.R) for a full example
using the bundled sample dataset (`inst/extdata/test.tab`, `test.gmt`).

## Functions

- `eem_search(expr, geneSets, ...)` — run the EEM search
- `read_gmt(file)` — read a GMT gene set file into a named list
- `module_activity(expr, result, normalize = TRUE, pvalueCutoff = 6)` — build
  a module x sample activity matrix from `eem_search()` results
