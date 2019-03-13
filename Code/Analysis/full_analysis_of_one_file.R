library(tools)
library(stringr)
library(bipartite)
library(igraph)
library(glue)
library(assertthat)
library(magrittr)
library(tidyverse)

source("strona_overlap_nestedness.R")

set.seed(0)

analyze_one_file <- function(dat_file, out_file, MAX_SIZE, OVERWRITE) {
  if (file.exists(out_file)) {
    already_calculated <-
      read_csv(out_file, col_types=str_c("ccc", str_c(rep("d", 33), collapse="")))
    successfully_calculated <- already_calculated %>%
      mutate_all(is.finite) %>%
      gather("metric", "value", everything()) %>%
      filter(value) %>%
      use_series(metric)
  } else {
    already_calculated <- tibble()
    successfully_calculated <- c()
  }
  # if this is an empirical file, it will have a weights column (that is
  # currently ignored). Randomized versions are unweighted by default
  edge_list <- read_csv(dat_file, col_names=FALSE,
                        col_types=ifelse(str_detect(dat_file, "/edgelists/"),
                                         "ccd", "cc"))
  g <- graph_from_data_frame(edge_list, directed=FALSE)
  V(g)$type <- V(g)$name %in% edge_list$X2
  B <- get.incidence(g)

  n_rows <- nrow(B)
  n_cols <- ncol(B)
  n_links <- sum(B)
  size <- n_rows * n_cols

  ## compute first few eigenvalues
  if (!all(c("l1", "l2", "l3", "ipr1", "ipr2", "ipr3") %in% successfully_calculated) | OVERWRITE) {
    spg <- tryCatch(igraph::spectrum(g, which=list(pos="LA", howmany=3), options=list(maxiter=5000)),
                    error=function(e) {return(list(values=NA, vectors=matrix(NA, 3, 3)))})
    if (is.na(spg$values[1]) & length(V(g)) < MAX_SIZE) {
      if (n_rows < n_cols) {
        spg$values <- eigen(B %*% t(B), only.values=TRUE)$values %>%
          sqrt() %>% sort(decreasing=TRUE) %>% head(3)
      } else {
        spg$values <- eigen(t(B) %*% B, only.values=TRUE)$values %>%
          sqrt() %>% sort(decreasing=TRUE) %>% head(3)
      }
    }
    l1 <- round(spg$values[1], 6)
    l2 <- round(spg$values[2], 6)
    l3 <- round(spg$values[3], 6)
    ipr1 <- sum(spg$vectors[,1]^4)
    ipr2 <- sum(spg$vectors[,2]^4)
    ipr3 <- sum(spg$vectors[,3]^4)
  } else {
    l1 <- already_calculated$l1
    l2 <- already_calculated$l2
    l3 <- already_calculated$l3
    ipr1 <- already_calculated$ipr1
    ipr2 <- already_calculated$ipr2
    ipr3 <- already_calculated$ipr3
  }

  deg <- as.vector(igraph::degree(g))
  d_row <- rowSums(B)
  d_col <- colSums(B)
  C_row <- (sum(d_row) - n_rows) / size
  C_col <- (sum(d_row) - n_cols) / size
  E_d_row <- 1 + C_row * (n_cols - 1)
  E_d_col <- 1 + C_col * (n_rows - 1)
  deg_het_row <- mean(d_row^2) / mean(d_row)^2
  deg_het_col <- mean(d_col^2) / mean(d_col)^2

  if (!all(c("l1_lap", "l2_lap", "l3_lap", "alg_conn") %in% successfully_calculated) | OVERWRITE) {
    lap_mat <- igraph::laplacian_matrix(g)
    lap_graph <- igraph::graph_from_adjacency_matrix(lap_mat, mode="undirected", weighted=TRUE)
    lap_mat_eigs <- tryCatch(spectrum(lap_graph, which=list(pos="LA", howmany=3),
                             options=list(maxiter=5000))$values,
                             error=function(e) {return(c(NA, NA, NA))})
    if (is.na(lap_mat_eigs[1]) & nrow(lap_mat) < MAX_SIZE) {
      lap_mat_eigs <- eigen(lap_mat, only.values=TRUE, symmetric=TRUE)$values}
    l1_lap <- lap_mat_eigs[1]
    l2_lap <- lap_mat_eigs[2]
    l3_lap <- lap_mat_eigs[3]
    alg_conn <- tryCatch(spectrum(lap_graph, which=list(pos="SA", howmany=2),
                         options=list(maxiter=5000))$values[2],
                         error=function(e) {return(NA)})
    if (is.na(alg_conn) & length(lap_mat_eigs) == nrow(lap_mat)) alg_conn <- tail(lap_mat_eigs, 2)[1]
  } else {
    l1_lap <- already_calculated$l1_lap
    l2_lap <- already_calculated$l2_lap
    l3_lap <- already_calculated$l3_lap
    alg_conn <- already_calculated$alg_conn
  }

  ## some motif metrics
  if (!all(c("H2", "H3", "H4", "H17") %in% successfully_calculated) | OVERWRITE) {
    # H1 <- length(E(g)) ## this is just the number of edges
    H3 <- ((d_row - 1) %*% B %*% (d_col - 1))[1]
    H2_r <- sum(choose(d_row, 2))
    H4_r <- sum(choose(d_row, 3))
    H17_r <- sum(choose(d_row, 4))
    H2_c <- sum(choose(d_col, 2))
    H4_c <- sum(choose(d_col, 3))
    H17_c <- sum(choose(d_col, 4))
    H2 <- H2_r + H2_c
    H4 <- H4_r + H4_c
    H17 <- H17_r + H17_c
  } else {
    H2 <- already_calculated$H2
    H3 <- already_calculated$H3
    H4 <- already_calculated$H4
    H17 <- already_calculated$H17
  }

  if (!all(c("clustering_c", "clustering_r") %in% successfully_calculated) | OVERWRITE) {
    if (H3 < (1000 * MAX_SIZE)) {
      clustering_c <- clustering_tm(as_edgelist(g))
      clustering_r <- clustering_tm(B %>% t() %>% graph_from_incidence_matrix() %>% as_edgelist())
    } else {
      clustering_c <- NA
      clustering_r <- NA
    }
  } else {
    clustering_c <- already_calculated$clustering_c
    clustering_r <- already_calculated$clustering_r
  }
  if (!("Q" %in% successfully_calculated) | OVERWRITE) {
    if (length(V(g)) < MAX_SIZE & min(dim(B)) > 5) {
      ## modularity
      Q <- computeModules(B)@likelihood
    } else {
      ## modularity
      Q <- NA
    }
  } else {
    Q <- already_calculated$Q
  }

  if (!all(c("N.nodf", "N.temp", "N.olap") %in% successfully_calculated) | OVERWRITE) {

    if (length(V(g)) < MAX_SIZE & min(dim(B)) > 5) {
      ## nestedness
      N.nodf <- nestednodf(B)$statistic["NODF"]
      N.temp <- nestedtemp(B)$statistic["temperature"]
      N.olap <- OverlapNestedness(B)
    } else {
      ## nestedness
      N.nodf <- NA
      N.temp <- NA
      N.olap <- NA
    }
  } else {
    N.nodf <- already_calculated$N.nodf
    N.temp <- already_calculated$N.temp
    N.olap <- already_calculated$N.olap
  }

  if (!all(c("deg_assort", "l1_er", "l1_cm", "l2_mp", "l3_mp", "l1_reg") %in%
          successfully_calculated) | OVERWRITE) {
    deg_assort <- igraph::assortativity.degree(g)
    d_row <- rowSums(B)
    d_col <- colSums(B)
    pr <- (n_links - n_rows) /  ((n_cols - 1) * n_rows)
    pc <- (n_links - n_cols) / ((n_rows - 1) * n_cols)
    l1_er <- sqrt((1 + (n_cols - 1) * pr * (3 + (n_cols - 2) * pr)) / (1 + (n_cols - 1) * pr) *
                    (1 + (n_rows - 1) * pc * (3 + (n_rows - 2) * pc)) / (1 + (n_rows - 1) * pc))
    ## expected spectral abscissa for configuration model random graph
    l1_cm <- sqrt((mean(d_col^2) / mean(d_col)) * (mean(d_row^2) / mean(d_row)))
    ## expected second, third eigenvalues if marczenko-pasteur
    C <- (sum(d_row) - l1_cm^2) / size
    l2_mp <- sqrt(n_rows * C * (1 + sqrt(n_cols/n_rows))^2)
    l3_mp <- l2_mp - sqrt(n_rows^(1/3) * pi^(2/3) * C * (1 - C) *
                            (1 + sqrt(n_cols/n_rows))^(4/3) /
                            (n_cols/n_rows)^(1/6))
    l1_reg <- mean(d_col) * n_cols / length(V(g)) + mean(d_row) * n_rows / length(V(g))
  } else {
    deg_assort <- already_calculated$deg_assort
    l1_er <- already_calculated$l1_er
    l1_cm <- already_calculated$l1_cm
    l2_mp <- already_calculated$l2_mp
    l3_mp <- already_calculated$l3_mp
    l1_reg <- already_calculated$l1_reg
  }

  ## centralities
  if (!all(c("cent_between", "cent_close", "cent_eigen", "diam", "mean_path_length") %in%
          successfully_calculated) | OVERWRITE) {
    cent_between <- centr_betw(g, directed=FALSE)$centralization
    cent_close <- centr_clo(g, mode="total")$centralization
    cent_eigen <- centr_eigen(g)$centralization
    diam <- diameter(g)
    mean_path_length <- mean_distance(g)
  } else {
    cent_between <- already_calculated$cent_between
    cent_close <- already_calculated$cent_close
    cent_eigen <- already_calculated$cent_eigen
    diam <- already_calculated$diam
    mean_path_length <- already_calculated$mean_path_length
  }

  ## get identifying information
  split_path <- str_split(out_file, "/")[[1]] %>% tail(3)
  if (!dir.exists(dirname(out_file))) dir.create(dirname(out_file), recursive=TRUE)

  ## add to results
  write_csv(tibble(
    randomization = split_path[1],
    name = file_path_sans_ext(split_path[3]),
    type = split_path[2],
    l1, l2, l3,
    ipr1, ipr2, ipr3,
    l1_cm, l1_er, l1_reg,
    l2_mp, l3_mp,
    l1_lap, l2_lap, l3_lap,
    alg_conn,
    clustering_c, clustering_r,
    deg_assort,
    Q, N.olap, N.temp, N.nodf,
    H2, H3, H4, H17,
    cent_between, cent_close, cent_eigen,
    diam, mean_path_length,
    deg_het_row, deg_het_col),
    path=out_file)
}

## If running this file directly, uncomment the lines below (starting with `# cArgs <- commandArgs(TRUE)`)
## and then call from terminal with the following arguments (in order):
##     {DATA FILE} -- an edgelist file (such as one of those provided in the /Data directory).
##     {OUTPUT FILE} -- the location of the output, usually similar to the data file, but in the /Results directory.
##     [MAX NETWORK SIZE] -- the maximum size of network to run the more computationally heavy
##                    analyses upon (this is implemented to limit time and memory load).
##                    This argument is optional and defaults to 100 (a very conservative value)
##                    if not provided.
##     [OVERWRITE] -- a flag indicating whether or not old results should be recalculated. This
##                    argument is optional and will default to FALSE if not provided.
## Examples:
## Rscript full_analysis_one_file.R ../../Data/edgelists/crimes/Chicago_crimes_2016-01-01.csv ../../Results/empirical/crimes/Chicago_crimes_2016-01-01.csv
## Rscript full_analysis_one_file.R ../../Data/edgelists/crimes/Chicago_crimes_2016-01-01.csv ../../Results/empirical/crimes/Chicago_crimes_2016-01-01.csv 1000 TRUE

# cArgs <- commandArgs(TRUE)
# dat_file <- cArgs[1]
# x <- assert_that(!is.na(dat_file), msg="The data file is a necessary argument")
# out_file <- cArgs[2]
# x <- assert_that(!is.na(out_file), msg="The output file is a necessary argument")
# MAX_SIZE <- cArgs[3]
# MAX_SIZE <- ifelse(is.na(MAX_SIZE), 100, as.integer(MAX_SIZE))
# x <- assert_that(!is.na(as.integer(MAX_SIZE)),
#                  msg=glue("unable to parse {cArgs[2]} into an integer matrix size limit"))
# OVERWRITE <- cArgs[4]
# OVERWRITE <- ifelse(is.na(OVERWRITE), FALSE, as.logical(OVERWRITE))
# x <- assert_that(!is.na(OVERWRITE),
#                  msg=glue("unable to parse {cArgs[3]} into a logical OVERWRITE value"))
# analyze_one_file(dat_file, out_file, MAX_SIZE, OVERWRITE)
