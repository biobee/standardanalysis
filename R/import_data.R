#' Import your data
#'
#' Give your biom file and mapping file to create a phyloseq object.
#' If the configuration file is not loaded, then a test object will be made.
#'
#'
#' @param name_biom_file the biom file, in the form of "config$biom"
#' @param name_mapping_file the mapping file, in the form of "config$mapping"
#'
#' @return a phyloseq object
#' @export
#'
#' @examples
import_data <- function(name_biom_file, name_mapping_file) {
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
	# Sorting on SampleID
	order <- sort(phyloseq::sample_names(biom))
	phyloseq::otu_table(biom) <- phyloseq::otu_table(biom)[,order]

	ps_obj <- phyloseq::merge_phyloseq(biom, mapping)

	phyloseq::tax_table(ps_obj) <- gsub("[A-Za-z]_[0-9]?_+","", phyloseq::tax_table(ps_obj))


	return(ps_obj)
}
