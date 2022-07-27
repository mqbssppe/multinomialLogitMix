



> package.skeleton(name = "multinomialLogitMixB", code_files = "multinomialLogitMix.R")

R> library(RcppArmadillo)
R> RcppArmadillo.package.skeleton(name = "multinomialLogitMix", code_files = "multinomialLogitMix.R", example_code=FALSE)

T> cp multinomialLogitMix.cpp multinomialLogitMix/src/multinomialLogitMix.cpp

meta edit description kai namespace

R> library(Rcpp)
R> compileAttributes("multinomialLogitMix")

meta edit to man

meta R CMD buil check etc




R CMD build multinomialLogitMix
R CMD check --as-cran multinomialLogitMix_1.0.tar.gz
R CMD INSTALL multinomialLogitMix


