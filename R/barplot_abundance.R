#' Make an abundance barplot
#'
#' @param pseq Your phyloseq object.
#' @param transf Transformation applied to otu counts. Default (and recommended) is relative (aka percentage / "compositional")
#' @param ... Other arguments to be passed on to both micrbiome::aggregate_rare AND to microbiome::plot_composition.
#'
#' @return An abundance barplot.
#' @export
#'
#' @examples barplot_abundance(ps, transf = "compositional", prevalence = 5, detection = 0.2, include.lowest = TRUE, sample.sort = NULL, otu.sort = NULL, average_by = NULL, group_by = NULL)

barplot_abundance <- function(pseq,
		transf = "compositional",
		...) {
  args <- list(...)

	ps_transf <- microbiome::transform(pseq, transf)
	#   merge all taxa that are detected rare
	# Arguments: level, detection, prevalence, include.lowest
	ps_rank <- microbiome::aggregate_rare(ps_transf, ...)

	plot_rank <- microbiome::plot_composition(ps_rank,
		plot.type = "barplot",
		verbose = FALSE, ...) +
		ggplot2::ylab("Relative abundance") +
	  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
			vjust = 0.85,
			hjust = 0)) +
	  ggplot2::scale_fill_brewer(args$level, palette = "Paired")

	if (!is.null(args$group_by)) {
	  plot_rank$data <-
	    base::merge(
	      plot_rank$data,
	      microbiome::meta(ps_rank),
	      sort = F,
	      all.x = T,
	      by.x = "Sample",
	      by.y = "X.SampleID"
	    )
	}

	return(plot_rank)
}
