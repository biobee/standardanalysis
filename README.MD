# Preliminary analysis pipeline

This package is a quick way to explore your data. It gives an overview of the raw data,
after which you can decide yourself on the sampling depth for subsampling (if you want to do that). It will provide alpha-diversity boxplots (Shannon, Simpson, etc.), relative abundance barplots, PCA, heatmaps, PERMANOVA and differential abundance analysis with Aldex2, ANCOM-BC2 and/or Maaslin2.
Use the option to adjust settings in the functions to your preference. See the microbiome and phyloseq packages (Microbiome: https://microbiome.github.io/tutorials/, Phyloseq: https://joey711.github.io/phyloseq/) for the (other) options. Else, just use the default settings. Feel free to change the content if you want to learn more / play around, you can find the code under the "standardanalysis/R" folder in the directory where your package is installed.

## Software pre-requisites:
- R - [Windows download](https://cran.r-project.org/bin/windows/base/)
- R - [Mac download](https://cran.r-project.org/bin/macosx/)
- R-studio - [Windows/Mac download](https://posit.co/download/rstudio-desktop/)


Dependencies (R-packages):
- These will all be handled inside of the Rmd file.

To install (with vignettes)
``` r
devtools::install_github("plasde/standardanalysis", build_vignettes = TRUE)
```

Input:
- BIOM file (.biom)
- Mapping file (.txt) (tab-delimited)


## Running the analysis

After installing R and R-studio, install the package from github with the command shown above (you may have to install and load the devtools package first).

You can either use the functions separately (functions can be seen in the vignette) or open the vignette itself and edit it to load your data. Open the vignette like this:
```{r}
vignette('standard-analysis-tutorial', package = 'standardanalysis')
```
Alternatively, open in web browser: https://github.com/plasde/standardanalysis/tree/main/doc/standard-analysis-tutorial.Rmd

If you have any problems, questions or things you'd like to see added just let me know by email d.w.k.van.der.plas@acta.nl (fastest) or the issue-board on github.
