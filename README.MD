# Preliminary analysis pipeline

This script is a quick way to explore your data. It gives an overview of the raw data,
after which you can decide yourself on the sampling depth. It will provide alpha-diversity boxplots (Shannon, Simpson), relative abundance barplots are shown and a PCA is performed.
Feel free to change any settings to your preference. See the microbiome and phyloseq packages (Microbiome: https://microbiome.github.io/tutorials/, Phyloseq: https://joey711.github.io/phyloseq/) for the options. Else, just use the default settings. If you do want to try different parameters but coding is too difficult for now, try out the R-Shiny App.

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
- Mapping file (.txt) (tab-delimited)
- BIOM file (.biom)
- config.yml (You can find it on this page, fill in some information about your data)
- scripts folder (feel free to change the content if you want to play around)


## Running the analysis

After installing R and R-studio, download the whole ACTA-Preventive-Dentistry directory.

![image_dl_zip](https://github.com/plasde/ACTA-Preventive-Dentistry/blob/main/assets/images/download_zip.png)

In the directory, under Preliminary-analysis you will find the prelim_analysis(.RMD) file which is the script. First, fill out (& save) the config(.yml) file according to your data and what you want to see. For the convenience I would recommend to COPY your data (biom file + mapping file) into the 'data' subdirectory that is in ACTA-Preventive-Dentistry.

With R and Rstudio installed, open the .rmd file, choose 'knit' and pick 'html'. R-studio will now make a html file with your analysis on it. Alternatively, you can run the code chunks seperately and/or change however you like.

![image_knitr](https://github.com/plasde/ACTA-Preventive-Dentistry/blob/main/assets/images/knit_R.png)

If you have any problems, questions or things you'd like to see added just let me know by email d.w.k.van.der.plas@acta.nl (fastest) or the issue-board on github.
