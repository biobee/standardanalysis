heatmap_make <- function(ps, detection, prevalence) {
	
	# First we make a selection of core ASVs/OTUs of our dataset 
	ps_core <- microbiome::core(ps,
		detection = detection,
		prevalence = prevalence,
		include.lowest = TRUE)

	# (Any) Log transform for decreasing range (aka bring values closer together for colours)
	# microbiome::transform gives pseudocount of "min val / 2"
  	ps_core_clr <- microbiome::transform(ps_core, transform = "clr")

	# Pretty heatmap making -
	# change fontsize + cellheight to what works for your amount of samples
  	pheatmap::pheatmap(phyloseq::otu_table(ps_core_clr),
  		main = "Heatmap, RCLR transformed",
		fontsize = 6,
		cellheight = 5)
}
