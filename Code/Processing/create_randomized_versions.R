library(igraph)
library(stringr)
library(tidyverse)
library(parallel)

source("../processing_data/from_incidence_matrix_to_biparite_graph.R")
source("curveball.R")

OVERWRITE <- FALSE
MAX_ATTEMPTS <- 100

## get the data files
data_files <- Sys.glob("../../edgelists/*/*.csv")

source("sparse_graph_starting_from_spanning_tree.R")

erdos_renyi_random_g <- function(g) {
    ## get B from the igraph object
    B <- get.incidence(g)
    ## take some measures to ensure connectedness:
    connected <- FALSE
    n_attempts <- 0
    while (!connected) {
        if (n_attempts > MAX_ATTEMPTS) {
            cat("\nReached", MAX_ATTEMPTS, "ER randomizations without finding a connected one")
            return(NULL)
        }
        n_attempts <- n_attempts + 1
        if (ecount(g) > max(dim(B)) * log(vcount(g))) {
            ## then the ER graph is almost surely connected
            ER_g <- sample_bipartite(nrow(B), ncol(B), m=ecount(g), type="gnm")
        } else {
            ## ensure no empty rows before assigning remaining links
            ER_g <- sparse_erdos_renyi(B)
        }
        ## double check that it is connected
        connected <- is.connected(ER_g)
    }
    cat("     ER")
    return(ER_g)
}
configuration_model_random_g <- function(g) {
    B <- get.incidence(g)
    # curveball works best if matrix is wider than long
    if (nrow(B) > ncol(B)) B <- t(B)
    connected <- FALSE
    n_attempts <- 0
    while (!connected) {
        CM_B <- curve_ball(B, 50)
        if (is.null(CM_B)) next()
        dimnames(CM_B) <- dimnames(B)
        CM_g <- graph_from_incidence_matrix(CM_B)
        ## double check that it is connected
        connected <- is.connected(CM_g)
        if (n_attempts > MAX_ATTEMPTS) {
            warning("\nReached ", MAX_ATTEMPTS, " CM randomizations without finding a connected one\n")
            return(NULL)
        }
        n_attempts <- n_attempts + 1
    }
    cat("     CM")
    if (nrow(B) > ncol(B)) {
        CM_B <- t(CM_B)
        return(graph_from_incidence_matrix(CM_B))
    } else {
        return(CM_g)
    }
}

## for each, produce two randomized versions: one Erdős–Rényi and one
## configuration model (each modified to ensure a single connected component)
## NOTE: don't use parallel (set mc.cores=1) if you want status updates
tmp <- mclapply(data_files, mc.cores=1, function(dat_file) {
    cat(dat_file)
    set.seed(0)
    edge_list <- read_csv(dat_file, col_names=FALSE, col_types="ccd")
    ## ignore the weights when randomizing
    g <- graph_from_data_frame(edge_list %>% select(X1, X2), directed=FALSE)
    V(g)$type <- V(g)$name %in% edge_list$X2
    original_g <- g
    ## randomize by scrambling links while ensuring a simple, fully connected graph
    filename_er <- str_replace(dat_file, "/edgelists/", "/edgelists_er/")
    ## ensure that the proper output directory exists
    dir.create(str_c(str_split(dirname(filename_er), "/")[[1]], collapse="/"),
               showWarnings=FALSE, recursive=TRUE)
    if (!file.exists(filename_er) | OVERWRITE) {
        g <- erdos_renyi_random_g(original_g)
        if (!is.null(g)) g %>% as_edgelist() %>% as_tibble() %>% write_csv(path=filename_er, col_names=FALSE)
    } else cat("     --")
    ## randomize while maintaining degree distributions and ensuring a simple, fully connected graph
    filename_cm <- str_replace(dat_file, "/edgelists/", "/edgelists_cm/")
    dir.create(str_c(str_split(dirname(filename_cm), "/")[[1]], collapse="/"),
               showWarnings=FALSE, recursive=TRUE)
    if (!file.exists(filename_cm) | OVERWRITE) {
        g <- configuration_model_random_g(original_g)
        if (!is.null(g)) g %>% as_edgelist() %>% as_tibble() %>% write_csv(path=filename_cm, col_names=FALSE)
    } else cat("     --")
    cat("\n")
})

## completion monitoring:
data_frame(full_names=list.files("../..", ".csv", recursive=TRUE)) %>%
    filter(str_detect(full_names, "edgelists")) %>%
    separate(full_names, c("rand", "type", "net"), "/") %>%
    count(rand, type) %>%
    spread(rand, n) %>% #summarise_at(vars(contains("edgelists")), sum)
    filter(edgelists != edgelists_cm | edgelists != edgelists_er | edgelists_cm != edgelists_er)
