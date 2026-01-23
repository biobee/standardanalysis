#' Principal Component Analysis
#'
#' @param pseq Your phyloseq object
#' @param transf What transformation to apply to the data. Default is CLR, with shift of 1.
#' @param colour_by What variable to colour the points by. Default no colouring.
#' @param ellipse Whether to include an ellips (t-distribution), TRUE or FALSE.
#' @param longi_lines Whether you want to include lines to connect repeated measures, TRUE or FALSE.
#' @param ... Additional parameters to be passed on to autoplot.
#' @param pseudocount What pseudocount you want to use in case of (any) log-transformation.
#'
#' @return A PCA plot
#' @export
#' @import ggfortify
#'
#' @examples make_PCA(ps, "clr", pseudocount = 0.5)

make_pca <- function(pseq,
                     transf = 'clr',
                     colour_by = NULL,
                     ellipse = FALSE,
                     longi_lines = NULL,
                     pseudocount = 1,
                     ...) {
  # Transformation - automatically does shift for log transformations

  if (transf %in% c("clr", "alr", "log","log10")) {
    pseq <- microbiome::transform(pseq, "shift", shift = pseudocount)
  }

  pseq <- microbiome::transform(pseq, transf)

  # PCA calculation
  pr_comp <-
    stats::prcomp(t(phyloseq::otu_table(pseq)),
                  center = T,
                  scale. = F)

  # Eigenvalue scaling - by Bernd
  pr_comp$x <- t(t(pr_comp$x) / pr_comp$sdev)

  # Plotting
  p <- ggfortify:::autoplot.prcomp(
    pr_comp,
    colour = colour_by,
    data = microbiome::meta(pseq),
    scale = FALSE,
    ...
  )

  # Extra settings
  if (!is.null(longi_lines)) {
    p <-
      p + ggplot2::geom_path(ggplot2::aes(group = base::get(longi_lines)),
                             color = "grey")
  }

  if (ellipse == TRUE) {
    p <- p + ggplot2::stat_ellipse(
      type = "t",
      alpha = 0.3,
      linewidth = 1,
      ggplot2::aes(colour = base::get(colour_by))
      )
  }

  return(p)
}
