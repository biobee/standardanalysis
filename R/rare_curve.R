#' Rarecurve for estimating rarefaction depth
#'
#' @param pseq Your phyloseq object.
#' @param ... Other arguments to be passed to vegan::rarecurve
#'
#' @return A rarecurve.
#' @export
#'
#' @examples rare_curve(ps)

rare_curve <- function(pseq,
                       ...) {

	abu_table <- t(microbiome::abundances(pseq))
	vegan::rarecurve(abu_table,
			step = 100,
			col = c("blue", "darkred", "darkgreen", "darkorange", "darkmagenta", "darksalmon"),
			cex = 0.7,
			lwd = 2,
			xlab = "Sample depth",
			ylab = "Observed",
			main = "Rarefaction curve with observed OTUs",
			...
			)

	graphics::grid(nx = NULL,
		ny = NULL,
		lty = 2,      # Grid line type
		col = "gray", # Grid line color
		lwd = 1)
}
