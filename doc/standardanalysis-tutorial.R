## ----include = FALSE--------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup, include = FALSE-------------------------------------------------------------------------------------------------------------------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))
base::RNGkind(kind = "L'Ecuyer-CMRG") # for parallel RNG


library(standardanalysis)


## ----message = FALSE--------------------------------------------------------------------------------------------------------------------------------------
input_directory <- "Path/to/your/data"
biom_file_name <- "test_biom_water.biom"
mapping_file_name <- "test_mapping_file.txt"
treatment_variable <- "Source"

# setwd(input_directory)


# ps <- import_data(biom_file_name, mapping_file_name, factors = "Source")
ps <- import_data()

# Function for cleaning up taxonomy
ps <- clean_taxonomy(ps)


read_num <- phyloseq::sample_sums(ps)
print(sort(read_num))
plot(sort(read_num), main = "Sample read depth")


rare_curve(ps, label = T)


## ----error=TRUE-------------------------------------------------------------------------------------------------------------------------------------------
try({
unwanted_samples <- c("MB.30.NP.0.MB.DNA048") # replace with your own sample names
ps <- phyloseq::subset_samples(ps, !X.SampleID %in% unwanted_samples)

# Check if the right ones are removed
microbiome::meta(ps)$X.SampleID
})


ps_rare <- rarefy_single(ps, depth = 0, mc = 1000)


alpha_boxplot(ps_rare, x_grp = treatment_variable, a_index = "observed", na.rm = T)
alpha_boxplot(ps_rare, x_grp = treatment_variable, a_index = "diversity_shannon") +
  ggplot2::ylab("Shannon diversity index")


ps_list <- rarefy_multiple(ps, sample.size = 2000, iter = 10)


alpha_df <- calculate_alpha_df(ps_list, measures = c("Shannon", "Simpson", "Chao1"))

ps_multi_ss <- phyloseq::subset_samples(ps, microbiome::meta(ps)$X.SampleID %in% unique(alpha_df$X.SampleID))
phyloseq::sample_data(ps_multi_ss) <- cbind(phyloseq::sample_data(ps_multi_ss), calculate_average_alpha_ps(alpha_df))


microbiome::boxplot_abundance(ps_multi_ss, x = treatment_variable, y = "median_Shannon")
microbiome::boxplot_abundance(ps_multi_ss, x = "X.SampleID", y = "median_Shannon")


# Note that you dont use average values here but "raw" alpha-div
mt_alpha <- multiple_test_alpha(
  alpha_df,
  pseq = ps_multi_ss,
  alpha_div = "Shannon",
  variable = treatment_variable,
  method = "wilcox.test"
)

acat(mt_alpha)
# median(x), mean(X), summary(x), min(x) (minimum)
# For example - using the above:
# median(mt_alpha)


adonis_results <- multiple_permanova(
  ps_list,
  distance = "bray",
  variable = "Readdepth + DNA_Conc",
  permutations = 999,
  # longit = "Source",
  ps_ref = NULL # Provide a phyloseq object
)

permanova_average(adonis_results)


# https://microbiome.github.io/tutorials/Composition.html
# Give taxonomic level, detection, prevalence.

# Optional to split plot by a variable: use "group_by = treatment_variable

# For grouping more than two variables, use outside of the function:
# + facet_wrap(Var1~Var2, scales = "free_x")
# or alternatively for nested-groups
# + facet_wrap(~Var1 + Var2, scales = "free_x")

barplot_abundance(
  ps,
  level = "Family",
  detection = 5 / 100,
  prevalence = 0 / nrow(microbiome::meta(ps))
) #+ facet_wrap(~My_grouping_variable + My_other_grouping_variable, scales = "free_x")


# Optional to also give a different shape to points based on a variable: shape = "My_other_variable"
make_pca(ps,
  colour_by = treatment_variable,
  ellipse = T,
  longi_lines = "Source"
)


make_heatmap(ps,
  detection = 7,
  prevalence = 0.4
)

# If this plot is too big to be shown properly, you can try copying the line with 'heatmap_make(...)' and paste it in the console-window down below this window.
# Do play around with the detection and prevalence, these defaults are just values that worked with the test data: detection = 7, prevalence = 0.4


## ----eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------
# ps_subset_sick65 <- phyloseq::subset_samples(ps, HealthStatus == "Sick" & age > 65)


ps_shift <- microbiome::transform(ps, transform = "shift", shift = 1)
ps_clr <- microbiome::transform(ps_shift, transform = "clr")


vegan::adonis2(
  as.formula(paste(
    "phyloseq::distance(phyloseq::otu_table(ps_clr), method = 'euclidean') ~",
    treatment_variable
  )),
  data = microbiome::meta(ps_clr),
  permutations = 9999
)


ps_shift <- microbiome::transform(ps, transform = "shift", shift = 1)
ps_clr <- microbiome::transform(ps_shift, transform = "clr")

longitudinal_variable <- ""

perm <- permute::how(nperm = 9999)
dat <- microbiome::meta(ps_clr)
permute::setBlocks(perm) <- with(dat, microbiome::meta(ps_clr)[[longitudinal_variable]])

vegan::adonis2(
  as.formula(paste(
    "phyloseq::distance(phyloseq::otu_table(ps_clr), method='euclidean') ~",
    treatment_variable
  )),
  data = dat,
  permutations = perm
)


# # Calculating instances and clr transform
# x <- ALDEx2::aldex.clr(as.data.frame(phyloseq::otu_table(ps)),
#   microbiome::meta(ps)[[treatment_variable]],
#   mc.samples = 128,
#   denom = "all",
#   verbose = TRUE
# )
# 
# # Testing
# x.tt <- ALDEx2::aldex.ttest(x,
#   paired.test = FALSE,
#   verbose = FALSE
# )
# x.effect <- ALDEx2::aldex.effect(x,
#   CI = TRUE,
#   verbose = FALSE,
#   paired.test = FALSE
# )
# 
# # Plotting
# x.all <- data.frame(x.tt, x.effect)
# par(mfrow = c(1, 2))
# ALDEx2::aldex.plot(x.all,
#   type = "MA", all.pch = 19, all.cex = 0.4,
#   called.col = "red", called.pch = 20, called.cex = 0.6,
#   thres.line.col = "darkgrey", thres.lwd = 1.5,
#   test = "welch", rare.col = "black",
#   rare = 0, rare.pch = 20, rare.cex = 0.2
# )
# ALDEx2::aldex.plot(x.all, type = "MW", test = "welch")
# 
# # Effect interval does not cross 0, has effect > 1 (or -1), has significant test with error correction
# x.all_sign <- subset(x.all, (x.all$effect.low > 0 & x.all$effect.high > 0) |
#   (x.all$effect.low < 0 & x.all$effect.high < 0) |
#   abs(x.all$effect) > 1 |
#   x.all$wi.eBH <= 0.05)
# 
# # I will print only the columns that are more interesting
# x.all_int <- x.all_sign[-c(1, 2, 5:9, 13)]
# 
# # Sorting on effect size
# print(x.all_int[order(abs(x.all_int$effect), decreasing = TRUE), ])


# params <- list(
#   # Model formula parameters
#   fixed_terms = "Timepoint + SubjectID", # Fixed effects formula
#   random_terms = "(Timepoint | SubjectID)", # Random effects formula
#   grouping_variable = NULL, # Grouping variable
# 
#   # General parameters
#   n_cores = 16, # Number of CPU cores or clusters
#   verbose = TRUE, # Verbose output
# 
#   # Preprocessing parameters
#   taxonomy_level = NULL, # Taxonomic level to be tested
#   prevalence_cutoff = 0.1, # Prevalence cut-off
#   library_size_cutoff = 1000, # Library size cut-off
# 
#   # Structural zeros parameters
#   structural_zero = TRUE, # Detect structural zeros
#   lower_bound = TRUE, # Classify a taxon as a structural zero using its asymptotic lower bound
# 
#   # Statistical parameters
#   p_adj_method = "holm", # P-value (multiple testing) adjustment method
#   alpha = 0.05, # Significance level
#   iter = 10, # Number of REML AND EM iterations
#   bootstrap = 10, # Number of bootstrap samples, should be 100+
# 
#   # Multi-group test parameters: set to FALSE if no grouping_variable!
#   global = TRUE, # Perform global test
#   pairwise = TRUE, # Perform pairwise tests
#   dunnet = TRUE, # Perform Dunnett's test
#   trend = TRUE, # Perform trend test
# 
#   # Pseudocount sensitivity analysis parameters
#   pseudo_sens = TRUE, # Sensitivity analysis
#   s0_perc = 0.05 # -th percentile of std. error values for each fixed effect (microarray sig.- SAM)
# )
# 
# out <- run_ancombc2(
#   ps = ps,
#   params = params
# )
# display_ancombc2_results(out, selected_results = c("global", "trend"))


# # For both the asv/otu table as well as meta data
# otu <- t(phyloseq::otu_table(ps))
# meta_data <- data.frame(microbiome::meta(ps))
# 
# # You can specifiy different GLMs/normalizations/transforms. These are according to  Nearing 2021 settings.
# fit_data <- Maaslin2::Maaslin2(
#   otu,
#   meta_data,
#   output = "My_maaslin2_analysis",
#   transform = "AST",
#   fixed_effects = treatment_variable,
#   normalization = "TSS",
# )
# 
# res_ms2 <- fit_data$results
# res_ms2_clean <- subset(res_ms2, res_ms2$qval <= 0.05)
# res_ms2_int <- res_ms2_clean[-c(2, 3, 5, 7, 9, 10)]
# print(res_ms2_int)


# sign_aldex <- rownames(x.all_int)
# sign_ancom <- res_anc_int[, 1]
# sign_maaslin2 <- res_ms2_int[, 1]
# 
# sign_3 <- intersect(sign_aldex, intersect(sign_ancom, sign_maaslin2))
# as.data.frame(tax_table(ps)[sign_3, -c(1:4)])

