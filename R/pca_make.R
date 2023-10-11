pca_make <- function(ps, transf = 'clr',
colour_by = NULL, ellipse = FALSE, longi_lines = FALSE, ...) {

	# Tranformation -- microbiome::transform with CLR gives pseudocount of "min value / 2" 
	phyloseq::otu_table(ps) <- microbiome::transform(phyloseq::otu_table(ps), transf)
	pr_comp <- stats::prcomp(t(phyloseq::otu_table(ps)))

	# Plotting
	p <- ggplot2::autoplot(pr_comp, colour = colour_by,
		size = 3, main = "Principle Component Analysis", data = microbiome::meta(ps), ...)

	# Extra settings
	if (longi_lines == TRUE) {
		p <- p + ggplot2::geom_line(aes(group = microbiome::meta(ps)$subjid),
			color = "grey")
	}

	if (ellipse == TRUE) {
		p <- p + ggplot2::stat_ellipse(type = "t",
		alpha = 0.3,
		size = 1,
		aes_string(colour = colour_by))
		#geom = "polygon" bij aes_string(fill=y)
	}
	# print(summary(pr_comp))
	return(p)
}