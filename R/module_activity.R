#' Compute a module x sample activity matrix
#'
#' For each expression module returned by \code{\link{eem_search}}, computes
#' its "activity" in each sample as the mean expression of its module genes
#' (the \code{moduleGenes} column), giving a matrix with one row per module
#' and one column per sample.
#'
#' @param expr the same numeric matrix (genes x samples, with gene rownames)
#'   that was passed to \code{\link{eem_search}}.
#' @param result a data.frame returned by \code{\link{eem_search}} (must have
#'   \code{id} and \code{moduleGenes} columns).
#' @param normalize standardize each gene (row of \code{expr}) to mean 0, sd 1
#'   before averaging (default \code{TRUE}), matching the row-normalization
#'   \code{eem_search()} applies internally before searching for modules, so
#'   that genes on different expression scales contribute comparably to the
#'   module's activity score. Set to \code{FALSE} to average raw expression
#'   values instead.
#' @param pvalueCutoff only include modules whose \code{pvalue} (-log10 scale,
#'   same convention as \code{pvalue1Cutoff} in \code{\link{eem_search}}) is
#'   at least this value. Default \code{6} (i.e. p < 1e-6). Set to \code{NULL}
#'   to include every module in \code{result}.
#'
#' @return a numeric matrix with one row per included module (named by
#'   \code{result$id}) and \code{ncol(expr)} columns (named by
#'   \code{colnames(expr)}).
#'
#' @examples
#' expr_file <- system.file("extdata", "test.tab", package = "eemr")
#' gmt_file  <- system.file("extdata", "test.gmt", package = "eemr")
#' expr <- as.matrix(read.table(expr_file, header = TRUE, row.names = 1,
#'                               sep = "\t", check.names = FALSE))
#' geneSets <- read_gmt(gmt_file)
#' res <- eem_search(expr, geneSets)
#' activity <- module_activity(expr, res)
#' dim(activity)
#' # only keep modules significant at p < 0.01 (-log10(0.01) = 2)
#' activity <- module_activity(expr, res, pvalueCutoff = -log10(0.01))
#' @export
module_activity <- function(expr, result, normalize = TRUE, pvalueCutoff = 6) {
	if (!is.matrix(expr) || !is.numeric(expr)) {
		stop("expr must be a numeric matrix")
	}
	if (is.null(rownames(expr))) {
		stop("expr must have rownames (gene ids)")
	}
	if (!all(c("id", "moduleGenes") %in% names(result))) {
		stop("result must be a data.frame with 'id' and 'moduleGenes' columns, as returned by eem_search()")
	}

	if (!is.null(pvalueCutoff)) {
		if (!("pvalue" %in% names(result))) {
			stop("result must have a 'pvalue' column to apply pvalueCutoff")
		}
		result <- result[result$pvalue >= pvalueCutoff, , drop = FALSE]
		if (nrow(result) == 0) {
			warning("no modules pass pvalueCutoff = ", pvalueCutoff)
		}
	}

	if (normalize) {
		expr <- t(scale(t(expr)))
	}

	moduleGeneList <- strsplit(result$moduleGenes, ";", fixed = TRUE)

	activity <- t(vapply(moduleGeneList, function(genes) {
		genes <- intersect(genes, rownames(expr))
		if (length(genes) == 0) {
			rep(NA_real_, ncol(expr))
		} else if (length(genes) == 1) {
			expr[genes, ]
		} else {
			colMeans(expr[genes, , drop = FALSE])
		}
	}, numeric(ncol(expr))))

	rownames(activity) <- result$id
	colnames(activity) <- colnames(expr)
	activity
}
