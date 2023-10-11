rarefy <- function(ps, depth = "NA") {
    if (depth == "NA" || depth == 0) {
        print("You chose no rarefaction.")
    
    } else if (is.numeric(depth)) {
        print(paste("You chose rarefaction on a depth of", depth))
        ps <- phyloseq::rarefy_even_depth(ps, rngseed = 3004, sample.size = depth, replace = F, verbose = F)
    
    } else {
        print("You did not provide a rarefaction depth.
        Consider the tables and graphs above and try again.
        If you do not want to rarefaction, choose 'NA'.")
    }
    
    # Remove samples with lower than 1500 reads
    ps <- phyloseq::prune_samples(sample_sums(ps) > 1499, ps)
    return(ps)
}