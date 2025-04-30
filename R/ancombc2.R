#' @title ANCOMBC2 Wrappers
#'
#' @description Functions for running ANCOM-BC2 analyses and displaying results.
#'
#' @section Functions:
#' - `run_ancombc2()`: Runs ANCOM-BC2 analysis on microbiome data.
#' - `display_results()`: Displays the results in an interactive table.
#'
#' @name ancombc2_wrappers
NULL
#' @export
params <- list(
  # Model formula parameters
  fixed_terms = "Timepoint + SubjectID", # Fixed effects formula
  random_terms = "(Timepoint | SubjectID)", # Random effects formula
  grouping_variable = NULL, # Grouping variable

  # General parameters
  n_cores = 16, # Number of CPU cores or clusters
  verbose = TRUE, # Verbose output

  # Preprocessing parameters
  taxonomy_level = NULL, # Taxonomic level to be tested
  prevalence_cutoff = 0.1, # Prevalence cut-off
  library_size_cutoff = 1000, # Library size cut-off

  # Structural zeros parameters
  structural_zero = TRUE, # Detect structural zeros
  lower_bound = TRUE, # Classify a taxon as a structural zero using its asymptotic lower bound

  # Statistical parameters
  p_adj_method = "holm", # P-value (multiple testing) adjustment method
  alpha = 0.05, # Significance level
  iter = 10, # Number of REML AND EM iterations
  bootstrap = 10, # Number of bootstrap samples, should be 100+

  # Multi-group test parameters: set to FALSE if no grouping_variable!
  global = TRUE, # Perform global test
  pairwise = TRUE, # Perform pairwise tests
  dunnet = TRUE, # Perform Dunnett's test
  trend = TRUE, # Perform trend test

  # Pseudocount sensitivity analysis parameters
  pseudo_sens = TRUE, # Sensitivity analysis
  s0_perc = 0.05 # -th percentile of std. error values for each fixed effect (microarray sig.- SAM)
)

#' Run ANCOM-BC2 Analysis
#'
#' This function executes the ANCOM-BC2 analysis on a given phyloseq object using specified parameters.
#' ANCOM-BC2 is a differential abundance analysis method for microbiome data.
#'
#' @param ps A `phyloseq` object containing the microbiome data.
#' @param params A named list of parameters for the ANCOM-BC2 analysis. See Details for required parameters.
#'
#' @return A list containing the ANCOM-BC2 analysis results:
#' \itemize{
#'   \item \code{res_global}: Results from global tests.
#'   \item \code{res_pair}: Results from pairwise comparisons.
#'   \item \code{dunn}: Results from Dunnett's test.
#'   \item \code{res_trend}: Results from trend analysis.
#'   \item \code{res}: Primary analysis results.
#' }
#'
#' @details
#' The `params` list should contain the following elements:
#' \describe{
#'   \item{taxonomy_level}{Taxonomic level to be tested (e.g., "Genus").}
#'   \item{prevalence_cutoff}{Prevalence cut-off for filtering taxa.}
#'   \item{library_size_cutoff}{Library size cut-off for filtering samples.}
#'   \item{structural_zero}{Logical, whether to detect structural zeros.}
#'   \item{lower_bound}{Logical, whether to classify a taxon as a structural zero using its asymptotic lower bound.}
#'   \item{fixed_terms}{Formula specifying fixed effects.}
#'   \item{random_terms}{Formula specifying random effects.}
#'   \item{grouping_variable}{Grouping variable for multi-group tests.}
#'   \item{p_adj_method}{Method for p-value adjustment (e.g., "holm").}
#'   \item{alpha}{Significance level.}
#'   \item{global}{Logical, whether to perform a global test.}
#'   \item{pairwise}{Logical, whether to perform pairwise tests.}
#'   \item{dunnet}{Logical, whether to perform Dunnett's test.}
#'   \item{trend}{Logical, whether to perform a trend test.}
#'   \item{pseudo_sens}{Logical, whether to perform pseudocount sensitivity analysis.}
#'   \item{s0_perc}{Percentile of standard error values for sensitivity analysis.}
#'   \item{n_cores}{Number of CPU cores to use.}
#'   \item{verbose}{Logical, whether to print verbose output.}
#'   \item{iter}{Number of iterations for REML and EM algorithms.}
#'   \item{bootstrap}{Number of bootstrap samples.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- run_ancombc2(ps, params)
#' }
#'
#' @export
run_ancombc2 <- function(ps, params) {
  out <- ANCOMBC::ancombc2(
    data = ps,
    tax_level = params$taxonomy_level,
    prv_cut = params$prevalence_cutoff,
    lib_cut = params$library_size_cutoff,
    struc_zero = params$structural_zero,
    neg_lb = params$lower_bound,
    fix_formula = params$fixed_terms,
    rand_formula = params$random_terms,
    group = params$grouping_variable,
    p_adj_method = params$p_adj_method,
    alpha = params$alpha,
    global = params$global,
    pairwise = params$pairwise,
    dunnet = params$dunnet,
    trend = params$trend,
    pseudo_sens = params$pseudo_sens,
    s0_perc = params$s0_perc,
    n_cl = params$n_cores,
    verbose = params$verbose,
    iter_control = list(tol = 1e-2, max_iter = params$iter, verbose = TRUE),
    em_control = list(tol = 1e-5, max_iter = params$iter),
    lme_control = lme4::lmerControl(),
    mdfdr_control = list(fwer_ctrl_method = "holm", B = params$bootstrap),
    trend_control = list(
      contrast = list(matrix(c(1, 0, -1, 1),
        nrow = 2, byrow = TRUE
      )),
      node = list(2), solver = "ECOS", B = params$bootstrap
    )
  )
  return(out)
}

#' Display ANCOM-BC2 Results
#'
#' This function presents the results of an ANCOM-BC2 analysis in an interactive table.
#'
#' @param out The output object from `run_ancombc2` (ANCOMBC::ancombc2 result).
#' @param analyses A character vector specifying which results to display.
#'   Possible values include "global", "pairwise", "dunnett", and "trend".
#'   If NULL, the primary analysis results are displayed.
#'
#' @return A `DT::datatable` object displaying the selected results.
#'
#' @examples
#' \dontrun{
#' display_results(result, analyses = c("global", "pairwise"))
#' }
#'
#' @export
display_ancombc2_results <- function(out, analyses = NULL, html = FALSE) {
  # Initialize list to store all DT objects
  dt_list <- list()

  # Helper functions remain the same
  get_result <- function(out, result_name) {
    switch(result_name,
      "global" = out$res_global,
      "pairwise" = out$res_pair,
      "dunnett" = out$dunn,
      "trend" = out$res_trend,
      NULL
    )
  }

  get_caption <- function(result_name) {
    switch(result_name,
      "global" = "ANCOM-BC2 Global Test",
      "pairwise" = "ANCOM-BC2 Pairwise Comparison",
      "dunnett" = "ANCOM-BC2 Dunnett's Test",
      "trend" = "ANCOM-BC2 Trend Analysis",
      "primary" = "ANCOM-BC2 Primary Analysis"
    )
  }

  # Process primary results
  numeric_cols <- sapply(out$res, is.numeric)
  out$res[numeric_cols] <- lapply(out$res[numeric_cols], function(x) round(x, 3))

  dt_list[["primary"]] <- ifelse(is.null(html), out$res, DT::datatable(out$res, caption = get_caption("primary")))

  # Process selected additional results
  if (!is.null(analyses)) {
    for (result in analyses) {
      temp_res <- get_result(out, result)
      if (!is.null(temp_res)) {
        numeric_cols <- sapply(temp_res, is.numeric)
        temp_res[numeric_cols] <- lapply(temp_res[numeric_cols], function(x) round(x, 3))
        dt_list[[result]] <- ifelse(is.null(html), temp_res, DT::datatable(temp_res, caption = get_caption(result)))
      }
    }
  }

  # Return all DT objects in a list
  return(dt_list)
}