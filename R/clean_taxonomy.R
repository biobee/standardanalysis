#' Cleaning of taxonomy.
#'
#' @param pseq Your phyloseq object.
#'
#' @return Your phyloseq object with cleaned taxonomy table.
#' @export
#'
#' @examples ps_clean <- clean_taxonomy(ps)

clean_taxonomy <- function(pseq) {

  phyloseq::tax_table(pseq) <- gsub("[A-Za-z]_[0-9]?_+", "", phyloseq::tax_table(pseq))

  # Get the old taxonomy table
  old_tax_table <- data.frame(phyloseq::tax_table(pseq))
  # For every row except the first, check if NA. If so, replace the lower ranks with "Unclassified [highest known taxa]"
  new_tax_table <- apply(old_tax_table,
                         1,
                         function(y) {
                           for(i in 2:7) {
                             if(grepl("uncultured", y[i], ignore.case = T)) { # if any of the values in row are uncultured:
                               y[i:7] <- replicate((8-i), paste("Uncultured", y[(i-1)], sep = " "))
                             }
                             if(base::is.na(y[i]) || y[i] == "" ||
                                grepl("metagenome|unclassified", y[i], ignore.case = T)) { # if any of the values in row are NA:
                               y[i:7] <- replicate((8-i), paste("Uncl", y[(i-1)], sep = " "))
                             }
                             if(grepl("unidentified", y[i], ignore.case = T)) { # if any of the values in row are unidentified:
                               y[i:7] <- replicate((8-i), paste("Unidentified", y[(i-1)], sep = " "))
                             }
                           }
                           return(y)
                         }
  )

  phyloseq::tax_table(pseq) <- t(new_tax_table)

  return(pseq)
}
