# Function for replacing NA
clean_taxonomy <- function(pseq) {
  
  # Get the old taxonomy table
  old_tax_table <- data.frame(phyloseq::tax_table(pseq))
  # For every row except the first, check if NA. If so, replace the lower ranks with "Unclassified [highest known taxa]"
  new_tax_table <- apply(old_tax_table,
                          1,
                          function(y) {
                            for(i in 2:7) {
                              if(is.na(y[i]) || y[i] == "") { # if any of the values in row are NA:
                                y[i:7] <- replicate((8-i), paste("Unclassified", y[(i-1)], sep = " "))
                              }
                            }
                            return(y)
                          }
                         )
  return(new.pseq <- phyloseq::phyloseq(phyloseq::otu_table(pseq), phyloseq::tax_table(t(new_tax_table)), phyloseq::sample_data(pseq)))
}