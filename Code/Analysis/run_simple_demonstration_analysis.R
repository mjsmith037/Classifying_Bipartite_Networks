## NOTE: This file takes about 15 minutes and 9 GB of memory to run, this memory
## requirement can be reduced by removing the largest data files from analysis
## using the MAX_SIZE parameter below.

library(tools)
library(stringr)
library(igraph)
library(magrittr)
library(tidyverse)

## ensures replicability
set.seed(0)

## max size for computationally intensive calculations
MAX_SIZE <- 1000

## This file assumes the data files are stored  in the Data/edgelists* folders of
## this repository
results <- lapply(Sys.glob("../../Data/edgelists*"), function(directory) {
  randomization <- str_replace_all(basename(directory),
                                   c("edgelists_er" = "erdos-renyi",
                                     "edgelists_cm" = "configuration",
                                     "edgelists"    = "empirical"))
  files <- list.files(directory, recursive=TRUE)
  nfiles <- length(files)
  cat("\n")
  lapply(1:length(files), function(ii) {
    # progress reporting
    cat("\r", randomization, ii, "/", nfiles, "\r")
    # if this is an empirical edgelist, it will have a weights column (that is
    # currently ignored). Randomized versions are unweighted by default
    edge_list <- read_csv(str_c(directory, "/", files[ii]), col_names=FALSE,
                          col_types=ifelse(randomization == "empirical", "ccd", "cc"),
                          progress=FALSE)
    # get both a (bipartite) igraph object and a corresponding incidence matrix
    g <- graph_from_data_frame(edge_list, directed=FALSE)
    V(g)$type <- V(g)$name %in% edge_list$X2
    B <- get.incidence(g)
    # compute first few eigenvalues using `spectrum` if possible ...
    spg <- tryCatch(spectrum(g, which=list(pos="LA", howmany=3), options=list(maxiter=5000)),
                    error=function(e) {return(list(values=NA, vectors=matrix(NA, 3, 3)))})
    # ... and manually if not, but care with large matrices and memory considerations
    if (is.na(spg$values[1]) & length(V(g)) < MAX_SIZE) {
      if (nrow(B) < ncol(B)) {
        spg <- sqrt(eigen(B %*% t(B)))
      } else {
        spg <- sqrt(eigen(t(B) %*% B))
      }
    }
    # round the values for storage and to reduce machine error
    l1 <- round(spg$values[1], 12)
    l2 <- round(spg$values[2], 12)
    l3 <- round(spg$values[3], 12)
    # what is the expectated spectral abscissa under a connected Erdos-Renyi random graph
    nr <- nrow(B)
    nc <- ncol(B)
    d_row <- rowSums(B)
    d_col <- colSums(B)
    pr <- (sum(B) - nr) /  ((nc - 1) * nr)
    pc <- (sum(B) - nc) / ((nr - 1) * nc)
    l1_er <- sqrt((1 + (nc - 1) * pr * (3 + (nc - 2) * pr)) / (1 + (nc - 1) * pr) *
                    (1 + (nr - 1) * pc * (3 + (nr - 2) * pc)) / (1 + (nr - 1) * pc))
    # what is the expectated spectral abscissa under a configuration model random graph
    l1_cm <- sqrt((mean(d_col^2) / mean(d_col)) * (mean(d_row^2) / mean(d_row)))
    # what is the expected  second, third eigenvalues if Marchenko–Pastur
    C <- (sum(d_row) - l1_cm^2) / (length(d_col) * length(d_row))
    l2_mp <- sqrt(length(d_row) * C * (1 + sqrt(length(d_col)/length(d_row)))^2)
    l3_mp <- l2_mp - sqrt(length(d_row)^(1/3) * pi^(2/3) * C * (1 - C) *
                            (1 + sqrt(length(d_col)/length(d_row)))^(4/3) /
                            (length(d_col)/length(d_row))^(1/6))
    # compile all of these results into a single data frame
    tibble(randomization,
           type = dirname(files[ii]),
           name = file_path_sans_ext(basename(files[ii])),
           l1, l2, l3,
           l1_cm, l1_er, l2_mp, l3_mp)
  }) %>% bind_rows()
}) %>% bind_rows()
cat("\n")

## write the output dataframe to file
write_csv(results, path="../../Results/SimpleDemo_results.csv")
