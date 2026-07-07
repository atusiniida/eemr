#' Coherence-based Expression Module (EEM) search
#'
#' For each gene set in \code{geneSets}, finds the largest coherently
#' expressed subset ("expression module") among its genes and assigns it a
#' p-value: an approximate, hypergeometric-based p-value (\code{pvalue1}) is
#' first used to screen candidates, and an accurate permutation-based
#' p-value (\code{pvalue2}) is then calculated for those that pass the
#' screen. \code{pvalue} is the one actually used (\code{pvalue2} when
#' available, otherwise \code{pvalue1}); all p-values are on the -log10
#' scale, so larger values are more significant.
#'
#' @param expr a numeric matrix of expression values with genes in rows and
#'   samples in columns. Must have rownames (gene ids); each row is
#'   internally standardized (mean 0, sd 1) before the search.
#' @param geneSets a named list of character vectors, e.g.
#'   \code{list(SET1 = c("gene1", "gene2"), SET2 = c("gene3", "gene4"))}.
#'   Names are used as gene set ids. See \code{\link{read_gmt}} to load these
#'   from a GMT file.
#' @param relativeRadius target fraction of all genes in \code{expr} that
#'   fall within the module search radius (default \code{0.05}).
#' @param pvalue1Cutoff cutoff (on -log10 scale) for the approximate p-value
#'   used to screen candidate gene sets before the accurate, permutation-based
#'   p-value is calculated. Default is \code{-log10(0.05)}. Use \code{-1} (or
#'   any negative value) to disable screening and calculate accurate p-values
#'   for every gene set.
#' @param itr number of permutations used to calculate the accurate p-value
#'   (default \code{300}).
#' @param minGeneSetSize,maxGeneSetSize gene sets whose size, after
#'   intersecting with the genes present in \code{expr}, falls outside
#'   \code{[minGeneSetSize, maxGeneSetSize]} are dropped (defaults \code{10}
#'   and \code{2000}).
#' @param coreGeneSize number of genes closest to the module center used to
#'   refine that center (default \code{10}).
#'
#' @return a data.frame with one row per gene set that passed the search,
#'   columns \code{id}, \code{pvalue}, \code{pvalue1}, \code{pvalue2},
#'   \code{nSeedGenes}, \code{nModuleGenes}, \code{absoluteRadius},
#'   \code{relativeRadius}, and \code{seedGenes}/\code{moduleGenes} (gene ids
#'   collapsed with \code{";"}), sorted by decreasing \code{pvalue}.
#'
#' @examples
#' expr_file <- system.file("extdata", "test.tab", package = "eemr")
#' gmt_file  <- system.file("extdata", "test.gmt", package = "eemr")
#' expr <- as.matrix(read.table(expr_file, header = TRUE, row.names = 1,
#'                               sep = "\t", check.names = FALSE))
#' geneSets <- read_gmt(gmt_file)
#' res <- eem_search(expr, geneSets)
#' head(res[, c("id", "pvalue", "nSeedGenes", "nModuleGenes")])
#' @export
eem_search <- function(expr, geneSets,
                        relativeRadius = 0.05,
                        pvalue1Cutoff = -log10(0.05),
                        itr = 300,
                        minGeneSetSize = 10,
                        maxGeneSetSize = 2000,
                        coreGeneSize = 10) {
	if (!is.matrix(expr) || !is.numeric(expr)) {
		stop("expr must be a numeric matrix")
	}
	if (is.null(rownames(expr))) {
		stop("expr must have rownames (gene ids)")
	}
	if (!is.list(geneSets) || is.null(names(geneSets)) || any(names(geneSets) == "")) {
		stop("geneSets must be a named list of character vectors")
	}
	geneSets <- lapply(geneSets, as.character)

	storage.mode(expr) <- "double"

	df <- eem_search_cpp(expr, geneSets,
	                      as.double(relativeRadius),
	                      as.double(pvalue1Cutoff),
	                      as.integer(itr),
	                      as.integer(minGeneSetSize),
	                      as.integer(maxGeneSetSize),
	                      as.integer(coreGeneSize))
	df[order(-df$pvalue), , drop = FALSE]
}

#' Read a GMT gene set file
#'
#' Reads a tab-delimited GMT file (as used by MSigDB/GSEA: one gene set per
#' line, columns are id, description, gene1, gene2, ...) into a named list
#' suitable for \code{\link{eem_search}}.
#'
#' @param file path to a GMT file.
#' @return a named list of character vectors (gene set id -> gene ids).
#' @examples
#' gmt_file <- system.file("extdata", "test.gmt", package = "eemr")
#' geneSets <- read_gmt(gmt_file)
#' length(geneSets)
#' geneSets[[1]]
#' @export
read_gmt <- function(file) {
	lines <- readLines(file)
	lines <- lines[nzchar(lines)]
	fields <- strsplit(lines, "\t")
	sets <- lapply(fields, function(x) x[-c(1, 2)])
	names(sets) <- vapply(fields, `[[`, character(1), 1)
	sets
}
