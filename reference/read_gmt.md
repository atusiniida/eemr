# Read a GMT gene set file

Reads a tab-delimited GMT file (as used by MSigDB/GSEA: one gene set per
line, columns are id, description, gene1, gene2, ...) into a named list
suitable for
[`eem_search`](https://atusiniida.github.io/eemr/reference/eem_search.md).

## Usage

``` r
read_gmt(file)
```

## Arguments

- file:

  path to a GMT file.

## Value

a named list of character vectors (gene set id -\> gene ids).
