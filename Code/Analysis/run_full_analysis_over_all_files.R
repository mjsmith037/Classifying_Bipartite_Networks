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

## As all of the following command line arguments are optional, this file can
## be run interactively as well as from the terminal. Simply change the defaults
## as you prefer and source this file to run locally. If running from the
## terminal, the arguments are (in order):
##     [OUTPUT DIRECTORY] -- the location of the output files
##     [MAX NETWORK SIZE] -- the maximum size of network to run the more computationally heavy
##                    analyses upon (this is implemented to limit time and memory load).
##     [OVERWRITE] -- a flag indicating whether or not old results should be recalculated.
## Note that the data directory is assumed to be ../../Data/s
## Examples:
## Rscript run_full_analysis_over_all_files.R
## Rscript run_full_analysis_over_all_files.R ../../Results/ 1000 TRUE

cArgs <- commandArgs(TRUE)
out_dir <- cArgs[1]
if (is.na(out_dir)) out_dir <- "../../Results/"

MAX_SIZE <- cArgs[2]
if (is.na(MAX_SIZE)) MAX_SIZE <- 100
assert_that(!is.na(as.integer(MAX_SIZE)),
            msg=glue("unable to parse {cArgs[2]} into an integer matrix size limit"))
OVERWRITE <- as.logical(cArgs[3])
if (is.na(OVERWRITE)) OVERWRITE <- FALSE
assert_that(!is.na(OVERWRITE),
            msg=glue("unable to parse {cArgs[3]} into a logical OVERWRITE value"))

## get all of the data files
data_files <- Sys.glob("../../Data/edgelists*/*/*.csv")
## and run the analysis script on each of them. You can adjust the number of cores
## used for this with the `mc.cores` argument below, but take care with memory usage
## for large files. 
x <- mclapply(data_files, mc.cores=1, function(dat_file) {
    out_file <- str_c(out_dir, str_extract(dat_file, "\\w+/\\w+/[\\w-_\\.]+\\.csv") %>%
                          str_replace_all(c("edgelists_er" = "erdos-renyi",
                                            "edgelists_cm" = "configuration",
                                            "edgelists"    = "empirical")))
    cat(dat_file, "\n")
    # system(str_c("Rscript analyze_one_file.R", dat_file, out_file, MAX_SIZE, OVERWRITE, sep=" "))
    analyze_one_file(dat_file, out_file, MAX_SIZE, OVERWRITE)
})

## completion monitoring:
full_results %>%
  group_by(randomization) %>%
  do(select(., -randomization, -type, -name) %>%
       summarise_all(. %>% is.na() %>% sum())) %>%
  do(select(., -randomization) %>%
       select_if(. %>% sum() %>% is_greater_than(0))) %>%
  ungroup() %>%
  mutate_all(replace_na, replace=0)
