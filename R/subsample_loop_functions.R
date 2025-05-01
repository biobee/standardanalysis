#' @title Subsample-loop
#'
#' @description
#' This suite of functions provides a workflow for performing and analyzing multiple rarefactions (random subsampling) of a \code{phyloseq} object. Multiple rarefactions help assess the variability and robustness of diversity and community composition metrics to subsampling.
#'
#' Most functions in this workflow (except \code{permanova_average}) support parallel execution using the \code{future} framework via \code{future.apply::future_lapply()}. This allows you to take advantage of multiple CPU cores or nodes to speed up computation, which is especially useful for large datasets or high-performance computing (HPC) environments.
#'
#' \strong{Getting started:}
#' \enumerate{
#'   \item Install the required package for parallelization:
#'     \itemize{
#'       \item \code{install.packages("future.apply")}
#'     }
#'   \item Before running these functions, set up your parallel backend using the \code{future} package. For example, to use all available cores on your machine:
#'     \itemize{
#'       \item \code{future::plan("multisession", workers = parallel::detectCores())}
#'     }
#'   \item For more information on available parallelization strategies and how to choose the best plan for your system (e.g., multicore, cluster, etc.), see the \href{https://future.futureverse.org/reference/plan.html}{future package documentation}.
#'   \item You can control the number of parallel workers via the \code{workers} argument (or \code{mc.cores} for some backends).
#'   \item For reproducible results in parallel, set the random number generator kind to L'Ecuyer-CMRG.
#'   \item Example setup:
#'     \itemize{
#'       \item \code{future::plan("multisession", workers = 4)}
#'       \item \code{base::RNGkind("L'Ecuyer-CMRG")}
#'     }
#' }
#'
#' The main steps in the workflow are:
#' \enumerate{
#'   \item \code{rarefy_multiple}: Create multiple rarefied versions of your phyloseq object.
#'   \item \code{calculate_alpha_df}: Calculate alpha-diversity metrics for each rarefied object.
#'   \item \code{calculate_average_alpha_ps}: Summarize alpha-diversity metrics across rarefactions.
#'   \item \code{multiple_test_alpha}: Test for group differences in alpha-diversity across rarefactions.
#'   \item \code{multiple_permanova}: Run PERMANOVA on each rarefied object to assess group differences in community composition.
#'   \item \code{permanova_average}: Aggregate and summarize PERMANOVA results across rarefactions.
#' }
#'
#' See the Details section for further information on each function.
#' @describeIn rarefy_multiple Creating multiple subsampled (rarefied) \code{phyloseq} objects
#' @param pseq (1) Your \code{phyloseq} object (or OTU-table-not recommended).
#' @param sample.size Numeric. Your preferred sample size for rarefying.
#' @param iter Numeric. The amount of times you would like to rarefy (amount of comparisons).
#' @param replace Logical. Perform subsampling with replacement or without (\code{TRUE} or \code{FALSE}).
#' @param seeds Vector of random seeds that you want to use. Possible to leave empty and have seed range to be 1:iterations.
#' @param ... Additional arguments to  be passed on to \code{phyloseq::rarefy_even_depth}
#' @return (1) A list of \code{phyloseq} objects with length of amount of \code{iter}
#' @details
#' (1) \code{rarefy_multiple} performs multiple rarefactions (random subsampling without or with replacement) on a \code{phyloseq} object. For each iteration, the function applies \code{\link[phyloseq]{rarefy_even_depth}} using a specified or automatically generated random seed. The result is a list of rarefied \code{phyloseq} objects, each representing a different random subsample at the specified depth. Rarefaction parameters (depth and replacement) are stored as attributes in the output list for reference.
#' @export
#' @examples pseq_list <- rarefy_multiple(pseq, sample.size = 8000, iter = 100)
# based on https://github.com/vmikk/metagMisc/blob/master/R/phyloseq_mult_raref.R
rarefy_multiple <- function(pseq,
                            sample.size = NULL,
                            iter = 10,
                            replace = FALSE,
                            seeds = NULL,
                            ...) {
  # Sanity check for the random number generator
  if (!is.null(seeds)) {
    if (length(seeds) != iter) {
      stop("Error: lenght of 'seeds' should be the same
            as the number of iterations.\n")
    }
    if (length(seeds) != length(unique(seeds))) {
      warning(
        "Warning: Provided seeds are not unique
            which leads to the identical results of random sampling.\n"
      )
    }
    if (!isTRUE(all(seeds == floor(seeds)))) {
      stop("Error: Seeds must only contain integer values.\n")
    }
  }

  # Define rarefication depth
  if (is.null(sample.size)) {
    stop("Error: You have to give a subsampling depth.\n")
  }

  # Prepare seed values
  if (is.null(seeds)) {
    seeds <- seq(1, iter)
  }
  # Do multiple rarefaction in parallel
  res_apply <- future.apply::future_lapply(seeds, function(ps) {
    phyloseq::rarefy_even_depth(
      pseq,
      sample.size = sample.size,
      rngseed = ps,
      replace = replace,
      verbose = FALSE, ...
    )
  }, future.seed = TRUE) # whether to generate seeds internally for reproducibility


  # Add rarefaction parameters as attributes to the resulting list
  base::attr(res_apply, which = "rarefaction_depth") <- sample.size
  base::attr(res_apply, which = "rarefaction_replacement") <- replace

  return(res_apply)
}

# ==============================================================================

#' @describeIn rarefy_multiple Estimating alpha-diversities
#' @param pseq_list (2) The list of \code{phyloseq} objects received from rarefy_multiple function.
#' @param measures Character. The alpha-diversity measure(s) you want to compute.
#' @return (2) A \code{data.frame} with the alpha-diversities requested for each \code{phyloseq} object.
#' @details
#' (2) \code{calculate_alpha_df} calculates alpha-diversity metrics for each rarefied \code{phyloseq} object in the input list. It applies `phyloseq::estimate_richness` to each object, optionally for a user-specified set of diversity measures. The results are combined into a single data frame, with each row corresponding to a sample from a particular rarefaction. This allows downstream analyses of the distribution of alpha-diversity across rarefactions.
#' @export
#' @examples alpha_df <- calculate_alpha_df(pseq_list, measures = c("Shannon", "Simpson", "Chao1"))
calculate_alpha_df <- function(
    pseq_list,
    measures = NULL) {
  alpha_df <- data.frame()

  # Estimate alpha diversity in parallel
  alphas <- future.apply::future_lapply(pseq_list, function(pseq) {
    phyloseq::estimate_richness(pseq, measures = measures)
  },
  future.seed = TRUE
  )

  # Stack (row-wise) each element of the list in a dataframe
  # Initialize an empty list to store the data frames
  alpha_list <- list()

  # Loop through each element in the 'alphas' list
  for (i in seq_along(alphas)) {
    # Get the current element
    current_df <- alphas[[i]]

    # Create a new column with row names
    current_df$X.SampleID <- rownames(current_df)

    # Append the modified data frame to the list
    alpha_list[[i]] <- current_df
  }

  # Combine all data frames in the list into a single data frame
  alpha_df <- do.call(rbind, alpha_list)
  row.names(alpha_df) <- NULL

  return(alpha_df)
}


# ==============================================================================
#' @describeIn rarefy_multiple Calculating of the averages of alpha-diversities
#' @param alpha_dataframe (3) The \code{data.frame} with alpha-diversities received from \code{calculate_alpha_df}.
#' @param alpha_div Character. The alpha-diversity measure(s) that you want to average. Default is \code{NULL} which gives all average for all measures in the dataframe.
#' @param averagef Character. The name of average function you want to use: median, mean, min (minimum) or max (maximum). Default is median (also recommended).
#' @return (3) A dataframe with the average alpha-diversities.
#' @details
#' (3) \code{calculate_average_alpha_ps} summarizes alpha-diversity values across multiple rarefactions for each sample. For each alpha-diversity metric, it computes the specified summary statistic (median, mean, min, or max) across all rarefied datasets. The output is a data frame where each row corresponds to a sample and each column to a summarized diversity measure. This facilitates comparison of central tendencies or ranges of diversity estimates across samples.
#' @export
#' @examples average_alphas <- calculate_average_alpha_ps(alpha_df, alpha_div = c("Shannon", "Simpson"), averagef = "min")
calculate_average_alpha_ps <- function(alpha_dataframe,
                                       alpha_div = NULL,
                                       averagef = "median") {
  # Catch wrong method argument
  allowed.methods <- c("median", "mean", "min", "max")
  if (!averagef %in% allowed.methods) {
    stop(
      "Non-supported average function specified. Allowed methods are one of: ",
      paste(allowed.methods, collapse = ", ")
    )
  }
  average.fun <- match.fun(averagef)

  # As we rely on named alpha_div, check if this exists in the dataframe
  alpha_divs <- names(alpha_dataframe)
  alpha_divs <- alpha_divs[alpha_divs != "X.SampleID"] # continue only with diversity-measures
  if (!all(alpha_div %in% names(alpha_dataframe))) {
    stop(
      "You call a measure not in your dataframe. You can select: ",
      paste(alpha_divs, collapse = ", "),
      "."
    )
  } else if (is.null(alpha_div)) {
    alpha_div <- alpha_divs
  }

  # Create a dataframe with the averages, grouped by X.SampleID
  avg_df <- aggregate(. ~ X.SampleID, data = alpha_dataframe[, c("X.SampleID", alpha_div)], FUN = average.fun)

  # Set column names based on averaging function used and alpha diversity used
  names(avg_df)[-1] <- paste(averagef, alpha_div, sep = ".")

  # Set SampleID as row names
  rownames(avg_df) <- avg_df$X.SampleID
  return(avg_df[-1])
}

# utils::globalVariables("X.SampleID")


# Alpha div test ===============================================================

#' @describeIn rarefy_multiple Testing of alpha-diversities between groups
#' @param alpha_dataframe (4) The \code{data.frame} with alpha-diversities received from \code{calculate_alpha_df}.
#' @param pseq Your original \code{phyloseq} object for the metadata.
#' @param alpha_div Character. The alpha-diversity measures that you want to test.
#' @param variable Character. The variable in the metadata with the groups that you want to test.
#' @param method Character. which test should be performed.
#' @param pair.by Character. In case of paired analysis, which is variable is the data paired on.
#' @return (4) A vector with the p-values for each \code{phyloseq} object.
#' @details
#' (4) \code{multiple_test_alpha} tests for differences in alpha-diversity between groups, separately for each rarefied dataset. Supported statistical tests include t-test, Wilcoxon rank-sum, Kruskal-Wallis, and Friedman test, with optional support for paired data. For each rarefied dataset, the appropriate test is applied to compare groups defined by a metadata variable. The output is a vector of p-values, one for each rarefaction, which can be aggregated or visualized to assess the robustness of group differences to subsampling variation.
#' @export
#' @examples alpha_pvals <- multiple_test_alpha(alpha_df, pseq, alpha_div = "Shannon", variable = "HealthStatus", method = "wilcoxon.test", pair.by = "SubjectID")
#' @examples median(alpha_pvals)
multiple_test_alpha <- function(alpha_dataframe,
                                pseq, # one pseq object with same metadata as in list
                                alpha_div,
                                variable,
                                method = "wilcox.test",
                                pair.by = NULL) {
  # Catch wrong method call
  allowed.methods <- c("t.test", "wilcox.test", "kruskal.test", "friedman.test")
  if (!(method %in% allowed.methods)) {
    stop(
      "Non-supported method specified. Allowed methods are one of: ",
      paste(allowed.methods, collapse = ", ")
    )
  }
  test.func <- match.fun(method)

  meta_df <- microbiome::meta(pseq)

  # Are all requested variables in metadata?
  if (!all(c(variable, pair.by) %in% names(meta_df))) {
    stop_str_meta <- "At least one of your variables is not in the metadata.\n"
    stop(stop_str_meta)
  }
  # Are all requested alpha diversities in alpha_data_frame?
  if (!all(alpha_div %in% names(alpha_dataframe))) {
    stop_str_alpha <- "At least one of your alpha_divs is not in the alpha data-frame.\n"
    stop(stop_str_alpha)
  }
  # Can't use t-test / wilcox with more than 2 groups
  if (length(unique(meta_df[[variable]])) > 2 &&
    method %in% c("t.test", "wilcox.test")) {
    stop("T.test or wilcoxon can only handle 2 groups.
         Use kruskal.test for non-paired and friedman.test for paired data.\n")
  }
  if (method == "kruskal.test") {
    if (!is.null(pair.by)) {
      stop("Kruskal-Wallis cannot test in a paired manner.\n")
    }
    test.func <- .test_kruskal
  } else if (method == "friedman.test") {
    test.func <- .test_friedman
  }

  paired <- FALSE
  if (!is.null(pair.by)) {
    ID <- pair.by
    if (length(pair.by) != 1) {
      stop("Provide only 1 pair.by!\n")
    }
    # Make sure all data are paired and ordered on ID variable
    SID_count <- table(microbiome::meta(pseq)[[X.SampleID]])
    if (max(SID_count) < 2) {
      stop("Your data are not paired on ", ID)
    }
    SID_all <- names(SID_count[SID_count == max(SID_count)])
    ss_meta <- subset(microbiome::meta(pseq), base::get(ID) %in% SID_all)

    # Refactor Subject/ pairby from 99 to 81
    ss_meta[[X.SampleID]] <- factor(ss_meta[[X.SampleID]])

    # Sort the dataframe on ID and grouping variable
    # Only complete pairs are left hereafter
    sort_by <- c(pair.by, variable)
    meta_df <- ss_meta[do.call(order, ss_meta[sort_by]), ]
    # Used prune_samples as subset_samples gave unexplained error with passing on sample_names argument
    pseq <-
      phyloseq::prune_samples(phyloseq::sample_names(pseq) %in% rownames(meta_df), pseq)

    paired <- TRUE
  }
  # Filter data from alpha_dataframe based on pseq because pairing might have removed unpaired samples;
  # We do not rely on nrow() but on samples
  # Note: this filtering HAS to be done after pairing samples in case of a paired test
  if (!all(phyloseq::sample_names(pseq) %in% unique(alpha_dataframe$X.SampleID))) {
    stop("Sample names in your alpha div dataframe and pseq object do not match.\n")
  } else if (length(unique(alpha_dataframe$X.SampleID)) != phyloseq::nsamples(pseq)) {
    warning(
      "Your alpha div dataframe and pseq object do NOT contain the same number of samples.\n",
      "Since all samples from pseq are in alpha_frame, we will select the pseq samples from alpha_dataframe.\n",
      "If you are doing paired-testing this is expected.\n"
    )
    alpha_dataframe <- subset(alpha_dataframe, X.SampleID %in% phyloseq::sample_names(pseq))
  }
  # Set the splitting variable as factor
  if (!is.factor(meta_df[[variable]])) {
    meta_df[[variable]] <- as.factor(meta_df[[variable]])
  }

  # Perform testing on pseq objects in parallel
  pvals <- future.apply::future_lapply(
    seq(1, nrow(alpha_dataframe),
      by = nrow(meta_df)
    ), function(pseq) {
      # Grab one pseq set, and do test on it
      values_single_ps <- alpha_dataframe[c(pseq:(pseq + nrow(meta_df) - 1)), alpha_div]

      # Always use meta[[variable]] here as meta has been sorted on the pairing
      # If a paired t-test or Wilcoxon signed rank test are used, subset samples based on unique variable values
      if (paired == TRUE & method %in% c("t.test", "wilcox.test")) {
        test.func(
          Pair(
            values_single_ps[meta_df[[variable]] == unique(meta_df[[variable]])[1]],
            values_single_ps[meta_df[[variable]] == unique(meta_df[[variable]])[2]]
          ) ~ 1,
          data = meta_df,
          variable = variable,
        )$p.value
        # Otherwise (independent samples, Friedman or Kruskal-Wallis) use this function call
      } else {
        test.func(
          values_single_ps ~ meta_df[[variable]],
          data = meta_df,
          variable = variable,
          values_single_ps = values_single_ps,
          ID = ID
        )$p.value
      }
    },
    future.seed = TRUE
  )
  pvals <- as.numeric(pvals)
  return(pvals)
}

.test_kruskal <- function(formula, ...) {
  stats::kruskal.test(formula)
}

.test_friedman <- function(formula, ...) {
  args <- list(...)
  lhs <- deparse(formula[[2]])
  stats::friedman.test(args[[lhs]], args$data[[args$variable]], args$data[[args$ID]])
}

# PERMANOVA ====================================================================

#' @describeIn rarefy_multiple Run PERMANOVA on Multiple Rarefied Phyloseq Objects
#' @param pseq_list (5) List or vector of rarefied \code{phyloseq} objects (from \code{rarefy_multiple}).
#' @param distance Character. Distance metric to use (e.g., \code{"bray"}, \code{"aitchison"}).
#' @param variable Character. Metadata variable(s) to test (e.g., \code{"Group"} or \code{"Var1 + Var2"}).
#' @param permutations Integer. Number of permutations for PERMANOVA.
#' @param pseudocount Numeric. Pseudocount to add to OTU table for Aitchison distance (default: 1).
#' @param ps_ref Optional. Reference \code{phyloseq} object for metadata subsetting.
#' @param longit Optional. Name of metadata variable for paired/longitudinal analysis.
#' @return (5) A list of adonis2 results, one per rarefied dataset.
#' @details
#' (4) \code{multiple_permanova} runs PERMANOVA (\code{\link[vegan]{adonis2}} ) on each rarefied \code{phyloseq} object to test for group differences in community composition. The specified distance metric (e.g., Bray-Curtis, Aitchison) is computed for each rarefied dataset, and the PERMANOVA is performed using the provided metadata variable(s). If a longitudinal variable is specified, permutation blocks are set up accordingly. The output is a list of adonis2 result objects, one per rarefied dataset, allowing assessment of the consistency of PERMANOVA results across rarefactions.
#' @export
#' @examples
#' permanova_results <- multiple_permanova(pseq_list, distance = "aitchison", variable = "Source", permutations = 9999, longit = "SubjectID")
multiple_permanova <- function(pseq_list,
                               distance,
                               variable,
                               permutations,
                               pseudocount = 1,
                               ps_ref = NULL,
                               longit = NULL) {
  adonis_result_list <- list()

  if (distance == "bray") {
    warning("Using Bray-Curtis distance.
    Data are not log2-transformed in this function.\n")
  }
  if (!is.null(ps_ref)) {
    if (is(ps_ref, "phyloseq")) {
      dat <- microbiome::meta(ps_ref)
    } else {
      stop("Error: 'ps_ref' is not a valid phyloseq object.\n")
    }
  } else {
    dat <- microbiome::meta(pseq_list[[1]])
  }
  # Set up blocks and overwrite "permutations" if longitudinal testing
  if (!is.null(longit)) {
    warning("Testing with strata or longitudinal design.\n")
    permutations <- permute::how(nperm = permutations)
    permute::setBlocks(permutations) <- with(dat, dat[[longit]])
  } else {
    warning("Testing without strata a.k.a.
    testing with cross-sectional design.\n")
  }

  # Perform multiple PERMANOVAs in parallel: fit adonis2 on every pseq in pseq_list
  adonis_result_list <- future.apply::future_lapply(pseq_list, function(pseq) {
    # Apply pseudo inside of loop because of different otu-tables per subsample
    if (distance == "aitchison") {
      phyloseq::otu_table(pseq) <- phyloseq::otu_table(pseq) + pseudocount
    }
    vegan::adonis2(
      as.formula(
        paste(
          "vegan::vegdist(t(phyloseq::otu_table(pseq)),",
          quote(distance), # vegdist requires quoted method arg
          ") ~ ",
          variable
        )
      ),
      data = dat,
      permutations = permutations
    )
  }, future.seed = TRUE)
  return(adonis_result_list)
}


# PERMANOVA - results aggregation ==============================================
#' @describeIn rarefy_multiple Aggregate PERMANOVA Results Across Multiple Rarefactions
#' @param results_adonis (5) A vector or list of PERMANOVA results as received from the \code{multiple_permanova}.
#' @param averagef Character. The averaging function to use. One of \code{"median"}, \code{"mean"}, \code{"min"}, or \code{"max"}.
#'   Default (and recommended) is \code{"median"}.
#' @param plot Logical. if \code{TRUE} (default), boxplots of pseudo-F statistics and p-values across rarefactions are displayed.
#' @return (5) A named list containing:
#'   1. \code{F.<averagef>}{The aggregated pseudo-F statistic (e.g., median pseudo-F).
#'   2. \code{F.IQR}The interquartile range (IQR) of pseudo-F statistics across rarefactions, rounded to 1 decimal, formatted as a string.
#'   3. \code{p.<averagef>} The aggregated p-value (e.g., median p-value).
#'   4. \code{p.Cauchy} The aggregated p-value computed using the Aggregated Cauchy Association Test (ACAT).
#'   5. \code{Fvals} Vector of pseudo-F statistics from each rarefaction.
#'   6. \code{pvals} Vector of p-values from each rarefaction.
#' }
#' @details
#' (5) \code{permanova_average} aggregates PERMANOVA results obtained from \code{\link[vegan]{adonis2}} on multiple rarefied \code{phyloseq} objects, the output of \code{multiple_permanova}.
#' It computes the summary statistic (median, mean, min, or max) of the PERMANOVA pseudo-F statistics and p-values, and combines p-values using the ACAT method (\code{\link{acat}}) for a robust overall significance measure. If \code{plot = TRUE} (default), boxplots visualize the distribution of pseudo-F statistics and p-values, with the ACAT p-value indicated by a red dashed line.
#'
#' @references
#' Liu, Y., & Xie, J. (2019).
#' Cauchy Combination Test: A Powerful Test With Analytic p-Value Calculation Under Arbitrary Dependency Structures.
#' \emph{Journal of the American Statistical Association}, \bold{115}(529), 393–402.
#' \doi{10.1080/01621459.2018.1554485}
#' @export
#' @examples
#' \dontrun{
#' results <- permanova_average(results_adonis, averagef = "median", plot = TRUE)
#' }
permanova_average <- function(
    results_adonis,
    averagef = "median",
    plot = TRUE) {
  allowed.methods <- c("median", "mean", "min", "max")

  if (!averagef %in% allowed.methods) {
    stop(
      "Non-supported average function specified. Allowed methods are one of: ",
      paste(allowed.methods, collapse = ", ")
    )
  }
  average.fun <- match.fun(averagef)
  # Initialise the vectors for F and p-values
  pvals <- c()
  fvals <- c()

  # Extract F and p-values from the adonis2 output list
  for (res in results_adonis) {
    pvals <- c(pvals, res$`Pr(>F)`[1]) # Extract p-values
    fvals <- c(fvals, res$F[1]) # Extract pseudo-F
  }
  # Calculate the F IQR (25th-75th percentile)
  F.IQR <- quantile(fvals, probs = c(0.25, 0.75))

  # Perform the Aggregated Cauchy Association Test (ACAT)
  p.Cauchy <- acat(pvals)

  # Create an output named list with the results
  out <- list(
    F = round(average.fun(fvals), 1),
    F.IQR = paste(round(F.IQR[1], 1), "-", round(F.IQR[2], 1)),
    p = average.fun(pvals),
    p.Cauchy = p.Cauchy,
    Fvals = fvals,
    pvals = pvals
  )

  names(out)[c(1, 3)] <- c(paste0("F.", averagef), paste0("p.", averagef))

  if (plot == TRUE) {
# ggplot2-based plotting
    library(ggplot2)

    # Pseudo-F boxplot
    df_f <- data.frame(statistic = fvals)
    p1 <- ggplot(df_f, aes(x = "", y = statistic)) +
      geom_boxplot(fill = "gray", width = 0.4) +
      labs(
        x = NULL, y = "pseudo-F statistic"
      ) +
      theme_classic() +
      theme(axis.text.x = element_blank())

    print(p1)

    # p-value boxplot with ACAT line
    df_p <- data.frame(p_value = pvals)
    p2 <- ggplot(df_p, aes(x = "", y = p_value)) +
      geom_boxplot(fill = "white", width = 0.4) +
      geom_hline(yintercept = p.Cauchy, color = "red", linetype = "dashed", size = 0.8) +
      labs(
        x = NULL, y = "p-value"
      ) +
      annotate(
        "text", x = 0, y = 0,
        label = paste("Cauchy p =", signif(p.Cauchy, 3)),
        color = "red", vjust = -0.5, hjust = -0.5
      ) +
      theme_classic() +
      theme(axis.text.x = element_blank())

    print(p2)
  }
  return(out)
}

# MISCELLANEOUS ================================================================
#' @title Aggregated Cauchy Association Test (ACAT)
#' @description
#' The Aggregated Cauchy Association Test (ACAT), also known as the Cauchy Combination Test (CCT), combines p-values that may be arbitrarily dependent into a single aggregated p-value.
#' The test statistic is defined as a weighted sum of Cauchy transformations of individual p-values.
#' Under arbitrary dependency structures, the tail of the null distribution of the test statistic can be approximated by a Cauchy distribution, enabling p-value calculation.
#' @param pvals Numeric vector of p-values to combine. All values must be strictly between 0 and 1.
#' @return A single aggregated p-value (numeric).
#' @references
#' Liu, Y., & Xie, J. (2019).
#' Cauchy Combination Test: A Powerful Test With Analytic p-Value Calculation Under Arbitrary Dependency Structures.
#' \emph{Journal of the American Statistical Association}, \bold{115}(529), 393–402.
#' \doi{10.1080/01621459.2018.1554485}
#' @export
#' @examples
#' combined_p <- acat(c(0.01, 0.03, 0.2))
acat <- function(pvals) {
  # Check: all p-values should be in (0, 1)
  if (any(pvals <= 0 | pvals >= 1)) {
    stop("All p-values must be strictly between 0 and 1.")
  }

  # Number of p-values
  K <- length(pvals)

  weights <- rep(1 / K, K)

  # Transform each p-value to Cauchy scale
  tans <- tan((0.5 - pvals) * pi)

  # Weighted sum
  T <- sum(weights * tans)

  # Compute combined p-value
  p.Cauchy <- 0.5 - atan(T) / pi
  return(p.Cauchy)
}