## See:
## Strona, Giovanni, et al. "A fast and unbiased procedure to randomize
##     ecological binary matrices with fixed row and column totals." Nature
##     communications 5 (2014): 4114.
## Carstens, Corrie J. "Proof of uniform sampling of binary matrices with fixed
##     row sums and column sums for the fast curveball algorithm." Physical
##     Review E 91.4 (2015): 042812.
## Carstens, Corrie Jacobien, Annabell Berger, and Giovanni Strona. "Curveball:
##     a new generation of sampling algorithms for graphs with fixed degree
##     sequence." arXiv preprint arXiv:1609.05137 (2016).

require(tidyverse)
require(magrittr)

my_sample <- function(x, size, replace=FALSE, prob=NULL) {
    if (missing(size)) size <- length(x)
    x[sample.int(length(x), size, replace, prob)]
}


## curve ball algorithm code slightly modified from:
## https://github.com/queenBNE/Curveball/blob/master/goodshuffle.R
## input: incidence matrix: an nxm matrix representing one of the symmetric off-
## diagonal blocks of an adjacency matrix of a bipartite graph
## and a multiplier for the number of swaps to conduct in randomizing
curve_ball <- function(incidence_matrix, rep_multiplier){
    n_rows <- dim(incidence_matrix)[1]
    n_cols <- dim(incidence_matrix)[2]
    # convert the incidence matrix to an adjacency list
    adjacency_list <- lapply(1:n_rows,
                             function(row) which(incidence_matrix[row,] == 1))
    
    for (rep in 1:(rep_multiplier * n_rows)){
        # pick two rows to participate in the trade
        participants <- my_sample(1:n_rows, 2)
        # and get each of their links
        a <- adjacency_list[[participants[1]]]
        b <- adjacency_list[[participants[2]]]
        # find which links they have in common
        shared_links <- intersect(a,b)
        # if a and b are not proper subsets
        if (length(shared_links) %>% is_in(c(length(a), length(b))) %>% not()) {
            # get the elements that are different between a and b (potential trades)
            diff_links <- setdiff(c(a,b), shared_links)
            # get the break point between a and b when concatenated
            L <- length(a) - length(shared_links)
            # initialize the replacement link vectors
            new_a <- a
            new_b <- b
            # only move on once a substantive trade has occurred
            while (setequal(new_a, a)) {
                # scramble the potential trades 
                diff_links <- my_sample(diff_links)
                # and distribute them between the new vectors
                new_a <- c(shared_links, diff_links[1:L])
                new_b <- c(shared_links, diff_links[(L+1):length(diff_links)])
            }
            adjacency_list[[participants[1]]] <- new_a
            adjacency_list[[participants[2]]] <- new_b
        }
    }
    randomized_matrix <- matrix(0, n_rows, n_cols)
    lapply(1:n_rows,
           function(row) randomized_matrix[row, adjacency_list[[row]]] <<- 1)
    return(randomized_matrix)
}