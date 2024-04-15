#' Importing data
#'
#' This function will import the data you've provided in the config file and make it into a phyloseq object.
#' In case no config file was provided, test data is loaded.
#'
#' @param name_biom_file The name of the biom file.
#' @param name_mapping_file The name of the mapping file.
#'
#' @return pseq A phyloseq object.
#' @export
#' @importFrom utils globalVariables
#'
#' @examples import_data(my_biom_file, my_mapping_file)

import_data <- function(name_biom_file, name_mapping_file) {


  # TODO // automatically load test data if none given
	if (!exists(quote(config))) {
		biom <- microbiome::read_phyloseq(
			system.file("extdata", "test_biom_file.biom", package = "standardanalysis"),
			type = "biom"
		)
		mapping <- phyloseq::import_qiime_sample_data(
			system.file("extdata", "test_mapping_file.txt", package = "standardanalysis")
		)

	} else {
		biom <- microbiome::read_phyloseq(paste0(config$input_directory, name_biom_file), type = "biom")
		mapping <- phyloseq::import_qiime_sample_data(paste0(config$input_directory, name_mapping_file))
	}

	pseq <- phyloseq::merge_phyloseq(biom, mapping)

	if (length(unique(c(nrow(mapping)),
	  phyloseq::nsamples(biom),
	  phyloseq::nsamples(pseq))) != 1) {
	  warning(paste("You have", nrow(mapping), "samples in mapping file... \n
	          and", phyloseq::nsamples(biom), "in your biom file...\n
	          and they combined to be", phyloseq::nsamples(pseq), "samples in your phyloseq object."))
	}

	return(pseq)
}

utils::globalVariables("config")
