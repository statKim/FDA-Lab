library(Rcpp)

# R source files
path <- "/Users/hyunsung/Desktop/Rproject/robfpca/R/"
flist <- list.files(path)
flist <- setdiff(flist, c("RcppExports.R","robfpca-package.R"))
for (fname in flist) {
  cat(paste0(fname, "\n"))
  source(paste0(path, fname))
}


# C++ source files
path <- "/Users/hyunsung/Desktop/Rproject/robfpca/src/"
flist <- list.files(path)
flist <- setdiff(flist, "RcppExports.cpp")
flist <- c("IRLS.cpp","cov_Mest.cpp")
for (fname in flist) {
  cat(paste0(fname, "\n"))
  sourceCpp(paste0(path, fname))
}
