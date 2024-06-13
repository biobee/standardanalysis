#' Making alpha-diversity boxplots
#'
#' @param pseq Your phyloseq object.
#' @param x_grp The variable of the groups you want to compare.
#' @param a_index The alpha-diversity you want to compare.
#' @param extra_var Optional: extra variables that you want to include in the comparison.
#' @param ... Other arguments passed on to microbiome::boxplot_alpha.
#'
#' @return A boxplot with alpha-diversities per group
#' @export
#'
#' @examples alpha_boxplot(ps, x_grp = "Source", "Shannon", extra_var = "Time_group")

alpha_boxplot <- function(pseq,
                          x_grp,
                          a_index,
                          extra_var = NULL, ...) {
  # na.rm = T removes the numerical values of alpha measure that are NA. Not the NA-values in variable group.
  p <- microbiome::boxplot_alpha(
    pseq,
    x_grp,
    a_index,
    ...
  )

  if (!is.null(extra_var)) {
    p <-
      p + ggplot2::facet_wrap(~base::get(extra_var))
  }

  return(p)
}
