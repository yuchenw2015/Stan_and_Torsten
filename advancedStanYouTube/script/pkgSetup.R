pkgs <- c('GGally', 'tidyverse',
          'bayesplot', 'loo', 
          'gdata', 'Hmisc', 
          "mvtnorm", 'sampling')

parentScriptDir <- getwd()  ## You may need to mod this to be in the top level of scriptDir
pkgDir <- file.path(parentScriptDir, "pkg")
dir.create(pkgDir)
libDir <- file.path(parentScriptDir, "lib")
dir.create(libDir)

.libPaths(libDir)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11 -DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION")

install.packages(c("devtools", "BH"),
                 lib=libDir,
                 contriburl=c(contrib.url("http://r-forge.r-project.org","source"),
                              contrib.url("https://cran.rstudio.com/","source")),
                 destdir=pkgDir,
                 type="source",
                 dependencies = c("Depends", "Imports", "LinkingTo"),
                 INSTALL_opts="--no-multiarch")
 
 library(devtools)
 devtools::install_github('metrumresearchgroup/TorstenHeaders')
 library(torstenHeaders)
 install_torsten()

install.packages(pkgs,
                 lib=libDir,
                 contriburl=c(contrib.url("http://r-forge.r-project.org","source"),
                              contrib.url("https://cran.rstudio.com/","source")),
                 destdir=pkgDir,
                 type="source",
                 dependencies = c("Depends", "Imports", "LinkingTo"),
                 INSTALL_opts="--no-multiarch")

