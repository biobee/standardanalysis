rare_curve <- function(ps, ...) {
	abutab <- t(microbiome::abundances(ps))
	vegan::rarecurve(abutab,
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
