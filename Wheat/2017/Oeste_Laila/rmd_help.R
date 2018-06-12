#===================================
# Preparation
#===================================
#--- set working directory ---#
setwd('~/Box Sync/ResearchProjects/DIFM')

library(rmarkdown)
# library(knitrBootstrap)

#===================================
# Render an Rmd file
#===================================
render('./Codes/Oeste_Laila/analysis_report.Rmd','html_document')

#===================================
# Run all the R code chunks in an rmd file
#===================================
runAllChunks <- function(rmd, envir=globalenv()){
  tempR <- tempfile(tmpdir = ".", fileext = ".R")
  on.exit(unlink(tempR))
  knitr::purl(rmd, output=tempR)
  sys.source(tempR, envir=envir)
}

runAllChunks('./Codes/Oeste_Laila/analysis_report.Rmd')
