library(tools)
library(bipartite)
library(igraph)
library(stringr)
library(tidyverse)

source("../Processing/sparse_graph_starting_from_spanning_tree.R")
source("../Processing/from_incidence_matrix_to_biparite_graph.R")
source("strona_overlap_nestedness.R")
source("../Processing/curveball.R")

set.seed(0)
cArgs <- commandArgs(TRUE)
full_name <- cArgs[1]
nrands <- cArgs[2]
if (is.na(nrands)) nrands <- 1000

MAX_ATTEMPTS <- 100

# results_file <- "../Results/full_results.csv"
# 
# #### read in the results file
# results_df <- read_csv(results_file,
#                        ## specify type for these columns because they vary between networks or cause errors in parsing
#                        col_types=cols(nrows = col_double(),
#                                       H2 = col_double(),
#                                       H3 = col_double(),
#                                       H4 = col_double(),
#                                       feature1 = col_character(),
#                                       feature2 = col_character(),
#                                       feature3 = col_character(),
#                                       feature4 = col_character())) %>%
#     ## clean up some names
#     mutate(type = tolower(type),
#            type = str_replace(type, "movies", "actor collaboration"),
#            type = str_replace(type, "ecologicalinteractions", "ecological interactions"))

# relative_nest_and_mod <- function(full_name, nrands=1000) {
    load(full_name)
    if (str_detect(full_name, "BipartiteMotifs")) B <- Data$B else B <- get.incidence(g)
    if (nrow(B) < ncol(B)) A <- t(B) %*% B else A <- B %*% t(B)
    ## empirical values
    empir <- data_frame(N.rho = eigen(A, symmetric=TRUE, only.values=TRUE)$values[1],
                        N.olap = OverlapNestedness(B),
                        N.temp = nestedtemp(B)$statistic["temperature"],
                        N.nodf = nestednodf(B)$statistic["NODF"],
                        Q = computeModules(B)@likelihood)
    ## mean of nrands randomizations (connected ER)
    means <- lapply(1:nrands, function(ii) {
        connected <- FALSE
        n_attempts <- 0
        while (!connected) {
            CM_B <- curve_ball(B, 25)
            dimnames(CM_B) <- dimnames(B)
            CM_g <- graph_from_incidence_matrix(CM_B)
            ## double check that it is connected
            connected <- is.connected(CM_g)
            if (n_attempts > MAX_ATTEMPTS) {
                cat("\nReached", MAX_ATTEMPTS, "CM randomizations without finding a connected one")
                return(NA)
            }
            n_attempts <- n_attempts + 1
        }
	B <- CM_B
#        B <- get.incidence(sparse_erdos_renyi(B))
        if (nrow(B) < ncol(B)) A <- t(B) %*% B else A <- B %*% t(B)
        data_frame(N.rho = eigen(A, symmetric=TRUE, only.values=TRUE)$values[1],
                   N.olap = OverlapNestedness(B),
                   N.temp = nestedtemp(B)$statistic["temperature"],
                   N.nodf = nestednodf(B)$statistic["NODF"],
                   Q = computeModules(B)@likelihood)
    }) %>% bind_rows()
    if (max(means %>% summarise_all(~sum(is.na(.)))) > 100) stop("less than 900 successful randomizations")
    means <- means %>% summarise_all(mean, na.rm=TRUE)
    # return((Empir-Mean)/Mean)
    write_csv((empir - means) / means, path=str_c("../Results/nest_mod_individual_results_cm/", basename(file_path_sans_ext(full_name)), ".csv"))
# }
# nest_mod_df <- results_df %>%
#     mutate(type = tolower(type)) %>%
#     filter(type == "mutualism" | type == "antagonism" | type == "ecologicalinteractions",
#            randomization == "None") %>%
#     rowwise() %>%
#     mutate(type = ifelse(type == "ecologicalinteractions", tolower(feature1), type)) %>%
#     group_by(name, type) %>%
#     do(., relative_nest_and_mod(.$name)) %>%
#     ungroup()
# write_csv(nest_mod_df, path="../Results/nest_mod_results.csv")
