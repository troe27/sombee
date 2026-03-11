#!/bin/bash
module load R/4.4.2-gfbf-2024a
mkdir $HOME/Rlibs
export R_LIBS_USER=$HOME/Rlibs


Rscript -e '.libPaths(c("~/Rlibs", .libPaths()));Sys.setenv(R_LIBS_USER = "~/Rlibs");Sys.setenv(TMPDIR = Sys.getenv("SCRATCH", unset = tempdir()));options(repos = c(CRAN = "https://cloud.r-project.org"));install.packages(c("remotes", "ggplot2","RColorBrewer","lsa", "lattice"));remotes::install_local("hdp", build_vignettes = FALSE, upgrade = "never")'
