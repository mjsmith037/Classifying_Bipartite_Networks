library(tools)
library(tidyverse)
library(parallel)
library(stringr)
library(bipartite)
library(igraph)

source("strona_overlap_nestedness.R")

set.seed(0)

## max size for computationally intensive calculations
MAX_SIZE <- 1000

# cArgs <- commandArgs(TRUE)
# directory <- cArgs[1]
# my_file <- cArgs[2]
Sys.time()
results <- lapply(c("../Data/data/", "../Data/data_er/", "../Data/data_cm/"), function(directory) {
    if (directory == "../Data/data/") randomization <- "None"
    if (directory == "../Data/data_er/") randomization <- "Erdos-Renyi"
    if (directory == "../Data/data_cm/") randomization <- "Configuration model"
    files <- list.files(directory, recursive = TRUE)
    nfiles <- length(files)
    cat("\n")
    ## care with using multiple cores, as large webs have high memory loads
    mclapply(1:length(files), mc.cores=1, function(i) {
        cat("\r", randomization, i, "/", nfiles, "\r")
        my_file <- files[i]
        load(str_c(directory, my_file))
        B <- get.incidence(g)
        
        ## compute first few eigenvalues
        spg <- tryCatch(spectrum(g, which=list(pos="LA", howmany=3), options=list(maxiter=5000)),
                        error=function(e) {return(list(values=NA, vectors=matrix(NA, 3, 3)))})
        if (is.na(spg$values[1]) & length(V(g)) < MAX_SIZE) {
            if (nrow(B) < ncol(B)) {
                spg <- sqrt(eigen(B %*% t(B)))
            } else {
                spg <- sqrt(eigen(t(B) %*% B))
            }
        }
        ## leading eigenvalues
        l1 <- round(spg$values[1], 6)
        l2 <- round(spg$values[2], 6)
        l3 <- round(spg$values[3], 6)

        nr <- nrow(B)
        nc <- ncol(B)
        d_row <- rowSums(B)
        d_col <- colSums(B)
        pr <- (sum(B) - nr) /  ((nc - 1) * nr)
        pc <- (sum(B) - nc) / ((nr - 1) * nc)
        l1_er <- sqrt((1 + (nc - 1) * pr * (3 + (nc - 2) * pr)) / (1 + (nc - 1) * pr) *
                          (1 + (nr - 1) * pc * (3 + (nr - 2) * pc)) / (1 + (nr - 1) * pc))
        ## expected spectral abscissa for configuration model random graph
        l1_cm <- sqrt((mean(d_col^2) / mean(d_col)) * (mean(d_row^2) / mean(d_row)))
        ## expected second, third eigenvalues if marczenko-pasteur
        C <- (sum(d_row) - l1_cm^2) / (length(d_col) * length(d_row))
        l2_mp <- sqrt(length(d_row) * C * (1 + sqrt(length(d_col)/length(d_row)))^2)
        l3_mp <- l2_mp - sqrt(length(d_row)^(1/3) * pi^(2/3) * C * (1 - C) *
                                  (1 + sqrt(length(d_col)/length(d_row)))^(4/3) /
                                  (length(d_col)/length(d_row))^(1/6))
        
        return(data_frame(
        # dir.create(str_c("../Results/individual_results/", basename(directory), "/", dirname(my_file)),
        #            showWarnings=FALSE, recursive=TRUE)
        # write_csv(data.frame(
            randomization = randomization,
            name = metadata$name,
            type = metadata$type,
            feature1 = as.character(metadata$feature_1),
            feature2 = as.character(metadata$feature_2),
            feature3 = as.character(metadata$feature_3),
            feature4 = as.character(metadata$feature_4),
            nrows = metadata$nrows,
            ncols = metadata$ncols,
            connectance = sum(B) / prod(dim(B)),
            nlinks = metadata$n_links,
            l1, l2, l3,
            l1_cm, l1_er, l2_mp, l3_mp))#,
            # path=str_c("../Results/individual_results/",
            #            basename(directory), "/",
            #            file_path_sans_ext(my_file), ".csv"))
    }) %>% bind_rows() %>% return()
}) %>% bind_rows()
cat("\n")
write_csv(results, path="../Results/SimpleDemo_results.csv")
Sys.time()
