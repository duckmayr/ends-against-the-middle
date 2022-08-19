## Code to install the exact version of the packages we used
install.packages("devtools", quiet = TRUE)
library(devtools, quietly = TRUE)
pkgs = list(
    c(pkg = "Rcpp",      version = "1.0.7"),
    c(pkg = "tibble",    version = "3.1.4"),
    c(pkg = "dplyr",     version = "1.0.7"),
    c(pkg = "tidyr",     version = "1.1.3"),
    c(pkg = "ggplot2",   version = "3.3.5"),
    c(pkg = "bggum",     version = "1.0.2"),
    c(pkg = "coda",      version = "0.19-4"),
    c(pkg = "MCMCpack",  version = "1.5-0"),
    c(pkg = "pROC",      version = "1.17.0.1"),
    c(pkg = "ltm",       version = "1.1-1"),
    c(pkg = "wnominate", version = "1.2.5"),
    c(pkg = "oc",        version = "1.2")

)
invisible(sapply(pkgs, function(x) {
    suppressMessages(install_version(
        package = x["pkg"],
        version = x["version"],
        repos = "http://cloud.r-project.org",
        quiet = TRUE,
        upgrade = "never"
    ))
}))
install("code/GGUM", quiet = TRUE, upgrade = "never")
