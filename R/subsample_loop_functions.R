#' @title Subsample-loop
#'
#' @description These functions are part of a method where multiple different subsamples are taken from a phyloseq object, to give an idea as to what the spread is subsampled sets.
#' Because these functions are highly dependent on each other they are put together in one description. Different functions are numbered in the description.
#'
#' @describeIn rarefy_multiple Making multiple subsampled phyloseq objects
#' @param pseq Your phyloseq object (or OTU-table (not recommended)).
#' @param sample.size Your preferred sample size for rarefying.
#' @param iter The amount of times you would like to rarefy (amount of comparisons).
#' @param replace Perform subsampling with replacement or without (TRUE or FALSE).
#' @param seeds Vector of random seeds that you want to use. Possible to leave empty and have seed range to be 1:iterations.
#' @param ... Additional arguments to  be passed on to phyloseq::rarefy_even_depth
#'
#' @return (1) A list of phyloseq objects with length of amount of iterations.
#' @export
#'
#' @examples ps_list <- rarefy_multiple(ps, sample.size = 8000, iter = 100)
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
    seeds <- 1:iter
  }

  res_apply <- lapply(
    seeds,
    function(z, ...) {
      phyloseq::rarefy_even_depth(
        pseq,
        rngseed = z,
        sample.size = sample.size,
        replace = replace,
        verbose = FALSE,
        ... # pass additional arguments to  phyloseq::rarefy_even_depth
      )
    }, ...
  )

  # Add rarefaction parameters as attributes to the resulting list
  base::attr(res_apply, which = "rarefaction_depth") <- sample.size
  base::attr(res_apply, which = "rarefaction_replacement") <- replace

  return(res_apply)
}

# ==============================================================================

#' @describeIn rarefy_multiple Calculating of alpha-diversities
#'
#' @param pseq_list (2) The list of phyloseq objects received from rarefy_multiple function.
#' @param measures The alpha-diversity measures that you want to compute.
#'
#' @return (2) A dataframe giving the alpha-diversities requested for each phyloseq object.
#' @export
#'
#' @examples alpha_df <- calculate_alpha_df(ps_list, measures = c("Shannon", "Simpson", "Chao1"))
calculate_alpha_df <- function(pseq_list, measures = NULL) {
  alpha_df <- data.frame()
  alphas <- furrr::future_map(pseq_list, ~ {
    df <- phyloseq::estimate_richness(.x, measures = measures)
  })
  alpha_df <- plyr::ldply(alphas, function(x) {
    ID <- rownames(x)
    cbind(x, ID)
  })
  return(alpha_df)
}


# ==============================================================================

#' @describeIn rarefy_multiple Calculating of the averages of alpha-diversities
#'
#' @param alpha_dataframe (3) The dataframe with alpha-diversities received from calculate_alpha_df function.
#' @param alpha_div The alpha-diversity measures that you want to average. Default is NULL which gives all average for all measures in the dataframe.
#' @param averagef The type of average function you want to use: median, mean, min (minimum) or max (maximum). Default is median (also recommended).
#'
#' @return (3) A dataframe with the average alpha-diversities.
#' @export
#'
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
  alpha_divs <- alpha_divs[alpha_divs != "ID"] # continue only with diversity-measures
  if (!all(alpha_div %in% names(alpha_dataframe))) {
    stop(
      "You call a measure not in your dataframe. You can select: ",
      paste(alpha_divs, collapse = ", "),
      "."
    )
  } else if (is.null(alpha_div)) {
    alpha_div <- alpha_divs
  }

  # alpha_div as vector of measures to loop over
  # avg_list <- list()
  avg_df <- dplyr::group_by(alpha_dataframe, ID) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(alpha_div), average.fun))
  names(avg_df)[-1] <- paste(averagef, alpha_div, sep = "_")
  avg_df <- as.data.frame(avg_df)
  return(avg_df)
}

# utils::globalVariables("ID")


# Alpha div test ===============================================================

#' @describeIn rarefy_multiple Testing of alpha-diversities between groups
#'
#' @param alpha_dataframe (4) The dataframe with alpha-diversities received from calculate_alpha_df function.
#' @param pseq Your original pseq object. This is required only for the metadata.
#' @param alpha_div The alpha-diversity measures that you want to test.
#' @param variable The variable in the metadata with the groups that you want to test.
#' @param method What test should be performed.
#' @param pair_by In case of paired analysis, what is variable is the data paired on.
#'
#' @return (4) A vector with the p-values for each phyloseq object.
#' @export
#'
#' @examples alpha_p_values <- multiple_test_alpha(alpha_df, pseq, alpha_div = "Shannon", variable = "HealthStatus", method = "wilcoxon.test", pair_by = "SubjectID")
#' @examples median(alpha_p_values)
multiple_test_alpha <- function(alpha_dataframe,
                                pseq,
                                # one pseq object with same metadata as in list
                                alpha_div,
                                variable,
                                method = "wilcox.test",
                                pair_by = NULL) {
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
  if (!all(c(variable, pair_by) %in% names(meta_df))) {
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
    if (!is.null(pair_by)) {
      stop("Kruskal-Wallis cannot test in a paired manner.\n")
    }
    test.func <- .test_kruskal
  } else if (method == "friedman.test") {
    test.func <- .test_friedman
  }

  paired <- FALSE
  if (!is.null(pair_by)) {
    ID <- pair_by
    if (length(pair_by) != 1) {
      stop("Provide only 1 pair_by!\n")
    }
    # Make sure all data are paired and ordered on ID variable
    SID_count <- table(meta_df[[ID]])
    if (max(SID_count) < 2) {
      stop("Your data are not paired on ", ID)
    }
    SID_all <- names(SID_count[SID_count == max(SID_count)])
    ss_meta <- subset(meta_df, base::get(ID) %in% SID_all)

    # Refactor Subject/ pairby from 99 to 81
    ss_meta[[ID]] <- factor(ss_meta[[ID]])

    # Sort the dataframe on ID and grouping variable
    # Only complete pairs are left hereafter
    sort_by <- c(pair_by, variable)
    meta_df <- ss_meta[do.call(order, ss_meta[sort_by]), ]
    # Used prune_samples as subset_samples gave unexplained error with passing on sample_names argument
    pseq <-
      phyloseq::prune_samples(phyloseq::sample_names(pseq) %in% rownames(meta_df), pseq)

    paired <- TRUE
  }

  # Filter data from alpha_dataframe based on pseq because pairing might have removed unpaired samples;
  # We do not rely on nrow() but on samples
  # Note: this filtering HAS to be done after pairing samples in case of a paired test
  if (!all(phyloseq::sample_names(pseq) %in% unique(alpha_dataframe$ID))) {
    stop("Sample names in your alpha div dataframe and pseq object do not match.\n")
  } else if (length(unique(alpha_dataframe$ID)) != phyloseq::nsamples(pseq)) {
    warning(
      "Your alpha div dataframe and pseq object do NOT contain the same number of samples.\n",
      "Since all samples from pseq are in alpha_frame, we will select the pseq samples from alpha_dataframe.\n",
      "If you are doing paired-testing this is expected.\n"
    )
    alpha_dataframe <- subset(alpha_dataframe, ID %in% phyloseq::sample_names(pseq))
  }
  # Set the splitting variable as factor
  if (!is.factor(meta_df[[variable]])) {
    meta_df[[variable]] <- as.factor(meta_df[[variable]])
  }

  list_of_test_p <- c()
  # For all samples in complete ps-list, with steps of X samples in one ps
  furrr::future_map(seq(1, nrow(alpha_dataframe), by = nrow(meta_df)), ~ {
    # Grab one ps set, and do test on it
    values_single_ps <- alpha_dataframe[c(.x:(.x + nrow(meta_df) - 1)), alpha_div]

    # Add p value to list of p-values - Treatment variable part of test can stay the same every iteration
    # Always use meta[[variable]] here as meta has been sorted on the pairing
    if (paired == T & method %in% c("t.test", "wilcox.test")) {
      list_of_test_p <- c(
        list_of_test_p,
        test.func(
          Pair(
            values_single_ps[meta_df[[variable]] == meta_df[[variable]][1]],
            values_single_ps[meta_df[[variable]] == meta_df[[variable]][2]]
          ) ~ 1,
          data = meta_df,
          variable = variable,
          values_single_ps = values_single_ps,
          ID = ID
        )$p.value
      )
    } else {
      list_of_test_p <- c(
        list_of_test_p,
        test.func(
          values_single_ps ~ meta_df[[variable]],
          data = meta_df,
          variable = variable,
          values_single_ps = values_single_ps,
          ID = ID
        )$p.value
      )
    }
  })
  return(list_of_test_p)
}

.test_kruskal <- function(formula, ...) {
  return(stats::kruskal.test(formula))
}

.test_friedman <- function(formula, ...) {
  args <- list(...)
  lhs <- deparse(formula[[2]])
  rhs <- deparse(formula[[3]])
  return(stats::friedman.test(args[[lhs]], args$data[[args$variable]], args$data[[args$ID]]))
}

# PERMANOVA ====================================================================

#' @describeIn rarefy_multiple Performing PERMANOVA
#'
#' @param list_of_ps (5) The list with subsampled phyloseq objects.
#' @param distance The distance metric to be used in the PERMANOVA.
#' @param variable The variable in the metadata with the groups that you want to test.In the case you want to test more than one variable use this "variable1 + variable2". Interactions can also be tested, to do so instead of "+" use "*".
#' @param permutations The amount of permutations that should be done by PERMANOVA, typically something like 999, 9999, 99999. (For reason behind this study PERMANOVA.)
#' @param pseudocount Pseudo count to be used in the case of Aitchison as distance metric. Default = 1.
#' @param longit In case of paired analysis, what is variable is the data paired on.
#'
#' @return (5) A vector with the outcomes from the permanova for each phyloseq object.
#' @export
#'
#' @examples permanova_results <- multiple_permanova(ps_list, distance = "aitchison", variable = "Source + Time", permutations = 9999, longit = "SubjectID")
multiple_permanova <- function(list_of_ps,
                               distance,
                               variable,
                               permutations,
                               pseudocount = 1,
                               longit = NULL) {
  adonis_result_list <- list()

  if (distance == "bray") {
    warning("Using Bray-Curtis. Data is not log2-transformed.\n")
  }

  # For metadata: take the first element in list_of_ps because meta data does not change
  dat <- microbiome::meta(list_of_ps[[1]])
  # Set up blocks and overwrite "permutations" if longitudinal testing
  if (!is.null(longit)) {
    warning("Testing with strata or longitudinal design.\n")
    permutations <- permute::how(nperm = permutations)
    permute::setBlocks(permutations) <- with(dat, dat[[longit]])
  } else {
    warning("Testing without strata a.k.a. testing with cross-sectional design.\n")
  }

  adonis_result_list <- furrr::future_map(list_of_ps, ~ {
    # Apply pseudo inside of loop because of different otu-tables per subsample
    if (distance == "aitchison") {
      phyloseq::otu_table(.x) <- phyloseq::otu_table(.x) + pseudocount
    }

    ado <- vegan::adonis2(
      as.formula(
        paste(
          "vegan::vegdist(t(phyloseq::otu_table(.x)),",
          quote(distance), # vegdist requires quoted method arg
          ") ~ ",
          variable
        )
      ),
      data = dat,
      permutations = permutations
    )
  })
  return(adonis_result_list)
}


# PERMANOVA - list of dataframes - picking p-values ============================

#' @describeIn rarefy_multiple Performing averaging on p-values from PERMANOVA
#'
#' @param results_adonis (6) The list with PERMANOVA results as received from the multiple_permanova function.
#' @param averagef The average function to be used. This can be "median", "mean", "min", "max". Default (and recommended) is median.
#'
#' @return (6) A dataframe with the average outcomes per variable from the PERMANOVAs.
#' @export
#'
#' @examples averages_PERMANOVA <- permanova_p_average(results_adonis, averagef = "median")
permanova_p_average <- function(results_adonis, averagef = "median") {
  allowed.methods <- c("median", "mean", "min", "max")
  if (!averagef %in% allowed.methods) {
    stop(
      "Non-supported average function specified. Allowed methods are one of: ",
      paste(allowed.methods, collapse = ", ")
    )
  }
  average.fun <- match.fun(averagef)

  p_values <- results_adonis |> furrr::future_map(~ {
    if (!is.na(.x$`Pr(>F)`[1])) {
      .x$`Pr(>F)`[1]
    }
  })
  # Coerce p_values to dataframe and transpose it
  res <- as.data.frame(p_values) |>
    t() |>
    as.data.frame()
  res[nrow(res) + 1, 1] <- average.fun(res[[1]])
  rownames(res) <- c(as.numeric(1:length(p_values)), paste(averagef))
  names(res) <- "P-values"
  return(res)
}
