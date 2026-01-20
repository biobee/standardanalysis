# Standard analysis R-package

This package is a quick way to explore your data. It gives an overview of the raw data,
after which you can decide yourself on the sampling depth for subsampling (if you want to do that). It will provide alpha-diversity boxplots (Shannon, Simpson, etc.), relative abundance barplots, PCA, heatmaps, PERMANOVA and differential abundance analysis with Aldex2, ANCOM-BC2 and/or Maaslin2.
Use the option to adjust settings in the functions to your preference. See the microbiome and phyloseq packages ([Microbiome](https://microbiome.github.io/tutorials/), [Phyloseq](https://joey711.github.io/phyloseq)) for the (other) options. Else, just use the default settings. Feel free to change the content if you want to learn more / play around, you can find the code under the "standardanalysis/R" folder in the directory where your package is installed.

## Software pre-requisites:
- R - [Windows download](https://cran.r-project.org/bin/windows/base/)
- R - [Mac download](https://cran.r-project.org/bin/macosx/)
- R-studio - [Windows/Mac download](https://posit.co/download/rstudio-desktop/)

Dependencies (R-packages):
- These will all be handled inside of the package.

## To install (with vignettes)
First make sure you have the package devtools installed and loaded. You will also need the BiocManager package to install both phyloseq and microbiome packages, as they are not available as a standard CRAN package.

``` r
install.packages('devtools')
library('devtools')

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Install Bioconductor dependencies
BiocManager::install(c("phyloseq", "microbiome"))
```

Now you can install the package from github. Choose build_vignettes = TRUE to be able to acces the vignettes.

``` r
devtools::install_github("biobee/standardanalysis", build_vignettes = TRUE)
library("standardanalysis")
```

In case, you get an error mentioning Rcmd.exe, make sure you have the [R directory added to your path.](https://www.architectryan.com/2018/03/17/add-to-the-path-on-windows-10/)
In my case, I need to add: "C:\Program Files\R\R-4.3.3\bin\x64"

## Running the analysis
Input:
- BIOM file (.biom)
- Mapping file (.txt) (tab-delimited)

After installing R and R-studio, install the package from github (shown above, install and load the devtools package first).

You can either use the functions separately (functions can be seen in the vignette) or open the vignette itself and edit it to load your data.
The structure of the vignette is how you would run your own analysis. Open the vignette like this:
```{r}
vignette('standardanalysis-tutorial', package = 'standardanalysis')
```
Alternatively, open in [web browser](https://html-preview.github.io/?url=https://github.com/biobee/standardanalysis/blob/main/doc/standardanalysis-tutorial.html).

A couple things to keep in mind:
- Consider what kind of normalisation you want for you data.
  - Generally, you use one kind of normalisation.
  - if your input data is already subsampled, there is no need to subsample again.
  - if you want to do multiple subsampling, start with raw data and do not perform CLR-transformation.
  - if you want to do CLR transformation, do not subsample but do remove samples having lower than 1000 reads.
- The points mentioned above also influence what distance metrics you use (Bray-Curtis, Euclidean, Aitchison).

If you have any problems, questions or things you'd like to see added just put it on the issue-board on github.
