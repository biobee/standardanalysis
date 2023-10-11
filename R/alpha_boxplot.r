alpha_boxplot <- function(ps, x_grp, a_index,
	more_var = FALSE, extra_var = NULL) {
	if (!exists(quote(config))) {
		x_grp = "Source"
	}

	p <- microbiome::boxplot_alpha(ps,
	x_grp,
	a_index,
	na.rm = TRUE,
	violin = FALSE,
	show.points = TRUE,
	fill.colors = NA)
	if (more_var == TRUE) {
		p <- p + ggplot2::facet_grid(cols = vars(microbiome::meta(ps)$extra_var))
	}
	return(p)
}