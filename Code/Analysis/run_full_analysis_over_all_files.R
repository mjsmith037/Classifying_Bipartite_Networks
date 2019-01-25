library(tools)
library(bipartite)
library(igraph)
library(parallel)
library(glue)
library(assertthat)
library(magrittr)
library(stringr)
library(tidyverse)

source("full_analysis_of_one_file.R")

cArgs <- commandArgs(TRUE)
out_dir <- cArgs[1]
if (is.na(out_dir)) out_dir <- "../../results/individual_results/"

MAX_SIZE <- cArgs[2]
if (is.na(MAX_SIZE)) MAX_SIZE <- 100
assert_that(!is.na(as.integer(MAX_SIZE)),
            msg=glue("unable to parse {cArgs[2]} into an integer matrix size limit"))
OVERWRITE <- as.logical(cArgs[3])
if (is.na(OVERWRITE)) OVERWRITE <- FALSE
assert_that(!is.na(OVERWRITE),
            msg=glue("unable to parse {cArgs[3]} into a logical OVERWRITE value"))

## get all of the data files
data_files <- Sys.glob("../../edgelists*/*/*.csv")
## and run the analysis script on each of them
x <- mclapply(data_files, mc.cores=4, function(dat_file) {
    out_file <- str_c(out_dir, str_extract(dat_file, "\\w+/\\w+/[\\w-_]+\\.csv") %>%
                          str_replace_all(c("edgelists_er" = "erdos-renyi",
                                            "edgelists_cm" = "configuration",
                                            "edgelists"    = "empirical")))
    cat(dat_file, "\n")
    # system(str_c("Rscript analyze_one_file.R", dat_file, out_file, MAX_SIZE, OVERWRITE, sep=" "))
    analyze_one_file(dat_file, out_file, MAX_SIZE, OVERWRITE)
})