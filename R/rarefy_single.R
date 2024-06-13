#' Rarefy otu counts
#'
#' @param pseq Your phyloseq object.
#' @param depth The depth you want to rarefy to.
#' @param mc The minimum count that a sample should have to not be thrown out.
#'
#' @return A phyloseq object where samples rarefied or samples below minimum count are thrown out.
#' @export
#'
#' @examples rarefy_single(ps, depth = 4000, mc = 1000)

rarefy_single <- function(pseq, depth = 0, mc = 1000) {
    if (depth == 0) {
        print("You chose not to rarefy. Samples with less than mc will be removed.")
        pseq <- phyloseq::prune_samples(phyloseq::sample_sums(pseq) >= mc, pseq)
    } else if (is.numeric(depth)) {
        print(paste("You chose rarefaction on a depth of", depth))
        pseq <- phyloseq::rarefy_even_depth(pseq, rngseed = 3004, sample.size = depth, replace = F, verbose = F)
    } else {
        print("You did not provide a rarefaction depth.
        Consider the sample counts and rarecurves and try again.
        If you do not want to rarefy, choose depth 0.")
    }

    return(pseq)
}
