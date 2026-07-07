# Coherence-based Expression Module (EEM) search

For each gene set in `geneSets`, finds the largest coherently expressed
subset ("expression module") among its genes and assigns it a p-value:
an approximate, hypergeometric-based p-value (`pvalue1`) is first used
to screen candidates, and an accurate permutation-based p-value
(`pvalue2`) is then calculated for those that pass the screen. `pvalue`
is the one actually used (`pvalue2` when available, otherwise
`pvalue1`); all p-values are on the -log10 scale, so larger values are
more significant.

## Usage

``` r
eem_search(
  expr,
  geneSets,
  relativeRadius = 0.05,
  pvalue1Cutoff = -log10(0.05),
  itr = 300,
  minGeneSetSize = 10,
  maxGeneSetSize = 2000,
  coreGeneSize = 10
)
```

## Arguments

- expr:

  a numeric matrix of expression values with genes in rows and samples
  in columns. Must have rownames (gene ids); each row is internally
  standardized (mean 0, sd 1) before the search.

- geneSets:

  a named list of character vectors, e.g.
  `list(SET1 = c("gene1", "gene2"), SET2 = c("gene3", "gene4"))`. Names
  are used as gene set ids. See
  [`read_gmt`](https://atusiniida.github.io/eemr/reference/read_gmt.md)
  to load these from a GMT file.

- relativeRadius:

  target fraction of all genes in `expr` that fall within the module
  search radius (default `0.05`).

- pvalue1Cutoff:

  cutoff (on -log10 scale) for the approximate p-value used to screen
  candidate gene sets before the accurate, permutation-based p-value is
  calculated. Default is `-log10(0.05)`. Use `-1` (or any negative
  value) to disable screening and calculate accurate p-values for every
  gene set.

- itr:

  number of permutations used to calculate the accurate p-value (default
  `300`).

- minGeneSetSize, maxGeneSetSize:

  gene sets whose size, after intersecting with the genes present in
  `expr`, falls outside `[minGeneSetSize, maxGeneSetSize]` are dropped
  (defaults `10` and `2000`).

- coreGeneSize:

  number of genes closest to the module center used to refine that
  center (default `10`).

## Value

a data.frame with one row per gene set that passed the search, columns
`id`, `pvalue`, `pvalue1`, `pvalue2`, `nSeedGenes`, `nModuleGenes`,
`absoluteRadius`, `relativeRadius`, and `seedGenes`/`moduleGenes` (gene
ids collapsed with `";"`), sorted by decreasing `pvalue`.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
expr <- matrix(rnorm(200 * 10), nrow = 200,
                dimnames = list(paste0("g", 1:200), paste0("s", 1:10)))
geneSets <- list(setA = paste0("g", 1:20), setB = paste0("g", sample(200, 20)))
eem_search(expr, geneSets)
} # }
```
