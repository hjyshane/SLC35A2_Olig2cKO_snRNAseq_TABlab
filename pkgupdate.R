sink("bioc_update.log", split = TRUE)


userlib <- "/igm/home/hxy008/R/x86_64-pc-linux-gnu-library/4.3"
.libPaths(c(userlib, setdiff(.libPaths(), userlib)))
repos <- BiocManager::repositories()  # Bioconductor + CRAN
old <- old.packages(lib.loc = userlib, repos = repos)

if (!is.null(old)) {
    options(Ncpus = parallel::detectCores())
    BiocManager::install(rownames(old), lib = userlib, ask = FALSE, update = TRUE)
}

sink()


# check
options(repos = BiocManager::repositories())
repos <- BiocManager::repositories()
lp <- .libPaths()

for (p in lp) {
    old <- old.packages(lib.loc = p, repos = repos)
    cat("\n===", p, "===\n")
    if (is.null(old)) {
        cat("All up to date here.\n")
    } else {
        cat(nrow(old), "packages out-of-date\n")
        print(head(rownames(old), 10))
    }
}
