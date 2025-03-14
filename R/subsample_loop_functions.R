#' @title Subsample-loop
#'
#' @description These functions are part of a method where multiple different subsamples are taken from a phyloseq object,
#' to give an idea as to what the spread is subsampled sets.
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
  }, future.seeds = TRUE) # whether to generate seeds internally for reproducibility


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

  # Estimate alpha diversity in parallel
  alphas <- future.apply::future_lapply(pseq_list, function(ps) {
    phyloseq::estimate_richness(ps, measures = measures)
  },
  future.seeds = TRUE
  )

  # Stack (row-wise) each element of the list in a dataframe
  # Initialize an empty list to store the data frames
  alpha_list <- list()

  # Loop through each element in the 'alphas' list
  for (i in seq_along(alphas)) {
    # Get the current element
    current_df <- alphas[[i]]

    # # Create a new column with row names
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
  names(avg_df)[-1] <- paste(averagef, alpha_div, sep = "_")

  # Set SampleID as row names
  rownames(avg_df) <- avg_df$X.SampleID
  return(avg_df[-1])
}

# utils::globalVariables("X.SampleID")


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
                                pseq, # one pseq object with same metadata as in list
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

  # Perform testing on ps objects in parallel
  list_of_test_p <- future.apply::future_lapply(seq(1, nrow(alpha_dataframe), by = nrow(meta_df)), function(ps) {
    # Grab one ps set, and do test on it
    values_single_ps <- alpha_dataframe[c(ps:(ps + nrow(meta_df) - 1)), alpha_div]

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
  }, future.seeds = TRUE)
  list_of_test_p <- as.numeric(list_of_test_p)
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
#' @param subset A string with the expression to use to subset the data
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
                               subset = NULL,
                               permutations,
                               pseudocount = 1,
                               longit = NULL) {
  adonis_result_list <- list()

  if (distance == "bray") {
    warning("Using Bray-Curtis distance. Data is not log2-transformed in this function.\n")
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

  # Perform multiple PERMANOVAs in parallel: fit adonis2 on every ps in ps_list
  adonis_result_list <- future.apply::future_lapply(list_of_ps, function(ps) {
    # In case of subsetting data based on conditions supplied as a string
    if (!is.null(subset)) {
      ps <- phyloseq::subset_samples(ps, eval(parse(text = subset)))
    }
    # Apply pseudo inside of loop because of different otu-tables per subsample
    if (distance == "aitchison") {
      phyloseq::otu_table(ps) <- phyloseq::otu_table(ps) + pseudocount
    }
    vegan::adonis2(
      as.formula(
        paste(
          "vegan::vegdist(t(phyloseq::otu_table(ps)),",
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

  # Extract the F-statistic p-values
  p_values <- results_adonis |> future.apply::future_lapply(function(p) {
    if (!is.na(p$`Pr(>F)`[1])) {
      p$`Pr(>F)`[1]
    }
  }, future.seeds = TRUE)
  # Coerce p_values to dataframe and transpose it
  res <- as.data.frame(p_values) |>
    t() |>
    as.data.frame()
  # Calculate the statistic requested (average.fun) and append it to the tail of the dataframe
  res[nrow(res) + 1, 1] <- average.fun(res[[1]])
  rownames(res) <- c(as.numeric(seq(nrow(res) - 1)), paste(averagef))
  names(res) <- "p-values"
  return(res)
}
