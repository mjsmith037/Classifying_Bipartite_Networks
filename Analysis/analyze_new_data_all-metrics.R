library(tools)
library(tidyverse)
library(parallel)
library(stringr)
library(bipartite)
library(igraph)

source("strona_overlap_nestedness.R")

set.seed(0)

## max size for computationally intensive calculations
MAX_SIZE <- 5000

cArgs <- commandArgs(TRUE)
directory <- cArgs[1]
my_file <- cArgs[2]

# results <- lapply(c("../Data/data/", "../Data/data_er/", "../Data/data_cm/"), function(directory) {
    if (directory == "../Data/data/") randomization <- "None"
    if (directory == "../Data/data_er/") randomization <- "Erdos-Renyi"
    if (directory == "../Data/data_cm/") randomization <- "Configuration model"
#     files <- list.files(directory, recursive = TRUE)
#     nfiles <- length(files)
#     cat("\n")
#     mclapply(1:nfiles, mc.cores=1, function(i) {
#         cat("\r", randomization, i, "/", nfiles, "\r")
#         my_file <- files[i]
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
        l1 <- round(spg$values[1], 6)
        l2 <- round(spg$values[2], 6)
        l3 <- round(spg$values[3], 6)
        ipr1 <- sum(spg$vectors[,1]^4)
        ipr2 <- sum(spg$vectors[,2]^4)
        ipr3 <- sum(spg$vectors[,3]^4)

        deg <- as.vector(degree(g))
        d_row <- rowSums(B)
        d_col <- colSums(B)
        C_row <- (sum(d_row) - length(d_row)) / (length(d_col) * length(d_row))
        C_col <- (sum(d_row) - length(d_col)) / (length(d_col) * length(d_row))
        E_d_row <- 1 + C_row * (length(d_col) - 1)
        E_d_col <- 1 + C_col * (length(d_row) - 1)
        deg_het_row <- mean(d_row^2) / mean(d_row)^2
        deg_het_col <- mean(d_col^2) / mean(d_col)^2
        lap_mat <- laplacian_matrix(g)
        lap_graph <- graph_from_adjacency_matrix(lap_mat, mode="undirected", weighted=TRUE)
        lap_mat_eigs <- tryCatch(spectrum(lap_graph, which=list(pos="LA", howmany=3), options=list(maxiter=5000)),
                                 error=function(e) {return(list(values=NA))})
        if (is.na(lap_mat_eigs[1]) & nrow(lap_mat) < MAX_SIZE) lap_mat_eigs <- eigen(lap_mat, only.values=TRUE)$values
        l1_lap <- lap_mat_eigs[1]
        l2_lap <- lap_mat_eigs[2]
        l3_lap <- lap_mat_eigs[3]
        alg_conn <- tryCatch(spectrum(lap_graph, which=list(pos="SA", howmany=2), options=list(maxiter=5000))$values[2],
                             error=function(e) {return(NA)})
        if (is.na(alg_conn) & length(lap_mat_eigs) == nrow(lap_mat)) alg_conn <- tail(lap_mat_eigs, 2)[1]
        ## some motif metrics
        # H1 <- length(E(g)) ## this is just the number of edges
        H3 <- (d_row - 1) %*% B %*% (d_col - 1)
        H2_r <- sum(choose(d_row, 2))
        H4_r <- sum(choose(d_row, 3))
        # H9_r <- 
        H17_r <- sum(choose(d_row, 4))
        # H18_r <- 
        H2_c <- sum(choose(d_col, 2))
        H4_c <- sum(choose(d_col, 3))
        # H9_c <- 
        H17_c <- sum(choose(d_col, 4))
        # H18_c <- 
        # C4 <- 
        # C6 <- 
        
        if (length(V(g)) < MAX_SIZE) {
            ## modularity
            Q <- computeModules(B)@likelihood
            ## nestedness
            N.nodf <- nestednodf(B)$statistic["NODF"]
            N.temp <- nestedtemp(B)$statistic["temperature"]
            N.olap <- OverlapNestedness(B)
        } else {
            ## modularity
            Q <- NA
            ## nestedness
            N.nodf <- NA
            N.temp <- NA
            N.olap <- NA
        }

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
        l1_reg <- mean(d_col) * length(d_col) / length(V(g)) + mean(d_row) * length(d_row) / length(V(g))
        
        ## centralities
        cent_between <- centr_betw(g, directed=FALSE)$centralization
        cent_close <- centr_clo(g, mode="total")$centralization
        cent_eigen <- centr_eigen(g)$centralization
        diam <- diameter(g)
        mean_path_length <- mean_distance(g)

        ## add to results
        # return(data.frame(
        dir.create(str_c("../Results/individual_results/", basename(directory), "/", dirname(my_file)),
                   showWarnings=FALSE, recursive=TRUE)
        write_csv(data.frame(
            randomization = randomization,
            name = metadata$name,
            type = metadata$type,
            feature1 = as.character(metadata$feature_1),
            feature2 = as.character(metadata$feature_2),
            feature3 = as.character(metadata$feature_3),
            feature4 = as.character(metadata$feature_4),
            nrows = metadata$nrows,
            ncols = metadata$ncols,
            connectance = metadata$connectance,
            nlinks = metadata$n_links,
            l1 = l1,
            l2 = l2, 
            l3 = l3,
            ipr1 = ipr1,
            ipr2 = ipr2,
            ipr3 = ipr3,
            l1_cm = l1_cm,
            l1_er = l1_er,
            l1_reg, l2_mp, l3_mp,
            l1_lap, l2_lap, l3_lap,
            alg_conn,
            Q, N.olap, N.temp, N.nodf,
            H2 = H2_c + H2_r, H3, H4 = H4_c + H4_r, H17 = H17_c + H17_r,
            cent_between, cent_close, cent_eigen, diam, mean_path_length,
            deg_het_row, deg_het_col,
            stringsAsFactors=FALSE),
            path=str_c("../Results/individual_results/",
                       basename(directory), "/",
                       file_path_sans_ext(my_file), ".csv"))
#    }) %>% bind_rows() %>% return()
# }) %>% bind_rows()
# cat("\n")

# write_csv(results, path="../Results/.csv")
