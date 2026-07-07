# Compute a module x sample activity matrix

For each expression module returned by
[`eem_search`](https://atusiniida.github.io/eemr/reference/eem_search.md),
computes its "activity" in each sample as the mean expression of its
module genes (the `moduleGenes` column), giving a matrix with one row
per module and one column per sample.

## Usage

``` r
module_activity(expr, result, normalize = TRUE, pvalueCutoff = 6)
```

## Arguments

- expr:

  the same numeric matrix (genes x samples, with gene rownames) that was
  passed to
  [`eem_search`](https://atusiniida.github.io/eemr/reference/eem_search.md).

- result:

  a data.frame returned by
  [`eem_search`](https://atusiniida.github.io/eemr/reference/eem_search.md)
  (must have `id` and `moduleGenes` columns).

- normalize:

  standardize each gene (row of `expr`) to mean 0, sd 1 before averaging
  (default `TRUE`), matching the row-normalization
  [`eem_search()`](https://atusiniida.github.io/eemr/reference/eem_search.md)
  applies internally before searching for modules, so that genes on
  different expression scales contribute comparably to the module's
  activity score. Set to `FALSE` to average raw expression values
  instead.

- pvalueCutoff:

  only include modules whose `pvalue` (-log10 scale, same convention as
  `pvalue1Cutoff` in
  [`eem_search`](https://atusiniida.github.io/eemr/reference/eem_search.md))
  is at least this value. Default `6` (i.e. p \< 1e-6). Set to `NULL` to
  include every module in `result`.

## Value

a numeric matrix with one row per included module (named by `result$id`)
and `ncol(expr)` columns (named by `colnames(expr)`).

## Examples

``` r
if (FALSE) { # \dontrun{
res <- eem_search(expr, geneSets)
activity <- module_activity(expr, res)
# only keep modules significant at p < 0.01 (-log10(0.01) = 2)
activity <- module_activity(expr, res, pvalueCutoff = -log10(0.01))
} # }
```
