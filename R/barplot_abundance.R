barplot_abundance <- function(ps_rared,
		transf = "compositional",
		rank = "Family",
		detection_rate = 5 / 100,
		prevalence_rate = 2/nrow(meta(ps_rared)), ...) {

	total_samples <- phyloseq::nsamples(ps_rared)
	ps_transf <- microbiome::transform(ps_rared, transf)
	#   merge all taxa that are detected rare
	ps_rank <- microbiome::aggregate_rare(ps_transf,
		level = rank,
		detection = detection_rate,
		prevalence = prevalence_rate)
		
	plot_rank <- microbiome::plot_composition(ps_rank,
		plot.type = "barplot",
		verbose = FALSE, ...) +
		ggplot2::xlab("Samples") +
		ggplot2::ylab("Relative abundance") +
		ggplot2::ggtitle("Abundance barplot") +
		ggplot2::theme(axis.text.x = element_text(angle = 90,
			vjust = 0.85,
			hjust = 0)) +
		ggplot2::scale_fill_brewer(palette = "Paired")

	return(plot_rank)
}
