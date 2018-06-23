library(igraph)
library(stringr)
library(tidyverse)
library(parallel)

source("../Processing/from_incidence_matrix_to_biparite_graph.R")
source("curveball.R")

OVERWRITE <- FALSE
MAX_ATTEMPTS <- 1000
set.seed(0)
## for randomizations that were unsuccessful above
# MAX_ATTEMPTS <- 10000
# set.seed(15987)

## get the data files
original_files <- list.files("../Data/data", recursive=TRUE, full.names=TRUE)
old_files <- list.files("../Data/old-data", recursive=TRUE, full.names=TRUE)

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
    connected <- FALSE
    n_attempts <- 0
    while (!connected) {
        # CM_B <- curve_ball(B, 25)
        CM_B <- sparse_configuration_model(B)
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
    return(CM_g)
}

## for each, produce two randomized versions: one Erdős–Rényi and one
## configuration model (each modified to ensure a single connected component)
## NOTE: don't use parallel (set mc.cores=1) if you want status updates
tmp <- mclapply(original_files, mc.cores=1, function(o_file) {
    cat(o_file)
    load(o_file)
    original_g <- g
    ## randomize by scrambling links while ensuring a simple, fully connected graph
    filename_er <- str_replace(o_file, "/data/", "/data_er/")
    ## ensure that the proper output directory exists
    dir.create(str_c(str_split(dirname(filename_er), "/")[[1]], collapse="/"),
               showWarnings=FALSE, recursive=TRUE)
    if (!file.exists(filename_er) | OVERWRITE) {
        g <- erdos_renyi_random_g(original_g)
        if (!is.null(g)) save(g, metadata, file=filename_er)
    } else cat("     --")
    ## randomize while maintaining degree distributions and ensuring a simple, fully connected graph
    filename_cm <- str_replace(o_file, "/data/", "/data_cm/")
    dir.create(str_c(str_split(dirname(filename_cm), "/")[[1]], collapse="/"),
               showWarnings=FALSE, recursive=TRUE)
    if (!file.exists(filename_cm) | OVERWRITE) {
        g <- configuration_model_random_g(original_g)
        if (!is.null(g)) save(g, metadata, file=filename_cm)
    } else cat("     --")
    cat("\n")
})

## for old-data
tmp <- mclapply(old_files, mc.cores=8, function(o_file) {
    cat(o_file)
    load(o_file)
    if (!str_detect(o_file, "igraph")) {
        original_g <- from_incidence_matrix_to_bipartite(Data$B)
        ## randomize by scrambling links while ensuring a simple, fully connected graph
        filename_er <- str_replace(o_file, "/old-data/", "/old-data_er/")
        ## ensure that the proper output directory exists
        dir.create(str_c(str_split(dirname(filename_er), "/")[[1]], collapse="/"),
                   showWarnings=FALSE, recursive=TRUE)
        if (!file.exists(filename_er) | OVERWRITE) {
            g <- erdos_renyi_random_g(original_g)
            Data <- list(B   = get.incidence(g),
                         RowDeg = Data$RowDeg,
                         ColDeg = Data$ColDeg,
                         NRows = Data$NRows,
                         NCols = Data$NCols,
                         NLinks = Data$NLinks,
                         Type = Data$Type,
                         SubType = Data$SubType,
                         BibTeX = Data$BibTeX,
                         Name = Data$Name,
                         RowIdentity = Data$RowIdentity,
                         ColIdentity = Data$ColIdentity)
            if (!is.null(g)) save(Data, file=filename_er)
        } else cat("     --")
        ## randomize while maintaining degree distributions and ensuring a simple, fully connected graph
        filename_cm <- str_replace(o_file, "/old-data/", "/old-data_cm/")
        dir.create(str_c(str_split(dirname(filename_cm), "/")[[1]], collapse="/"),
                   showWarnings=FALSE, recursive=TRUE)
        if (!file.exists(filename_cm) | OVERWRITE) {
            g <- configuration_model_random_g(original_g)
            Data <- list(B   = get.incidence(g),
                         RowDeg = Data$RowDeg,
                         ColDeg = Data$ColDeg,
                         NRows = Data$NRows,
                         NCols = Data$NCols,
                         NLinks = Data$NLinks,
                         Type = Data$Type,
                         SubType = Data$SubType,
                         BibTeX = Data$BibTeX,
                         Name = Data$Name,
                         RowIdentity = Data$RowIdentity,
                         ColIdentity = Data$ColIdentity)
            if (!is.null(g)) save(Data, file=filename_cm)
        } else cat("     --")
    } else {
        original_g <- g
        ## randomize by scrambling links while ensuring a simple, fully connected graph
        filename_er <- str_replace(o_file, "/old-data/", "/old-data_er/")
        ## ensure that the proper output directory exists
        dir.create(str_c(str_split(dirname(filename_er), "/")[[1]], collapse="/"),
                   showWarnings=FALSE, recursive=TRUE)
        if (!file.exists(filename_er) | OVERWRITE) {
            g <- erdos_renyi_random_g(original_g)
            if (!is.null(g)) save(g, metadata, file=filename_er)
        } else cat("     --")
        ## randomize while maintaining degree distributions and ensuring a simple, fully connected graph
        filename_cm <- str_replace(o_file, "/old-data/", "/old-data_cm/")
        dir.create(str_c(str_split(dirname(filename_cm), "/")[[1]], collapse="/"),
                   showWarnings=FALSE, recursive=TRUE)
        if (!file.exists(filename_cm) | OVERWRITE) {
            g <- configuration_model_random_g(original_g)
            if (!is.null(g)) save(g, metadata, file=filename_cm)
        } else cat("     --")

    }
    cat("\n")
})

## completion monitoring:
data_frame(full_names=list.files('..', '.RData', recursive=TRUE)) %>%
    filter(!str_detect(full_names, "rawdata")) %>%
    separate(full_names, c("rand", "type", "net"), "/") %>%
    mutate(rand = str_replace(rand, "old-", "")) %>%
    count(rand, type) %>%
    spread(rand, n) %>% #summarise_at(vars(contains("data")), sum)
    filter(data != data_cm | data != data_er | data_cm != data_er)
