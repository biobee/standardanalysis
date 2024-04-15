#' Make a heatmap
#'
#' @param pseq Your phyloseq object.
#' @param transf What transformation you want to use. Note that a log transformation is highly recommended.
#' @param pseudocount What pseudocount to use before log-transformation.
#' @param ... Other arguments to pass to microbiome::core. Such as prevalence and detection.
#'
#' @return A heatmap plot
#' @export
#' @import pheatmap
#'
#' @examples make_heatmap(ps, transf = "clr", pseudocount = 1, prevalence = 5, detection = 0.2, include.lowest = TRUE)

make_heatmap <-
  function(pseq,
           transf = "clr",
           pseudocount = 1,
           ...) {

    # First we make a selection of core ASVs/OTUs of our dataset
    ps_core <- microbiome::core(
      pseq,
      ...
    )

    # (Any) Log transform for decreasing range (aka bring values closer together for colours)
    if (transf %in% c("clr", "alr", "log","log10")) {
      ps_core <- microbiome::transform(ps_core, "shift", shift = pseudocount)
    }

    ps_core <- microbiome::transform(ps_core, transf)


    # Pretty heatmap making -
    # change fontsize + cellheight to what works for your amount of samples
    pheatmap::pheatmap(
      phyloseq::otu_table(ps_core),
      fontsize = 6,
      cellheight = 5
    )
  }
