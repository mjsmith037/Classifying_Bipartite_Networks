## a more consistent sampling function
my.sample <- function(x, size=1) x[sample(length(x), size)]

# create a random walk through the complete graph
# add a link every time a new node is discovered
# this will provide a skeleton for the graph: 
# a (random) spanning tree will connect all of the nodes
random_spanning_tree <- function(num_rows, num_cols) {
    tree_B <- matrix(0, num_rows, num_cols)
    discovered_rows <- rep(0, num_rows)
    discovered_cols <- rep(0, num_cols)
    my_row <- sample(1:num_rows, 1)
    my_col <- sample(1:num_cols, 1)
    discovered_rows[my_row] <- 1
    discovered_cols[my_col] <- 1
    tree_B[my_row, my_col] <- 1
    while(sum(c(discovered_cols, discovered_rows) == 0) > 0){
        my_row <- sample(1:num_rows, 1)
        if (discovered_rows[my_row] == 0){
            discovered_rows[my_row] <- 1
            tree_B[my_row, my_col] <- 1 
        }
        my_col <- sample(1:num_cols, 1)
        if (discovered_cols[my_col] == 0){
            discovered_cols[my_col] <- 1
            tree_B[my_row, my_col] <- 1 
        }
    }
    return(tree_B)
}

sparse_erdos_renyi <- function(B) {
  ## create a random spanning tree to ensure connectedness
  ER_B <- random_spanning_tree(nrow(B), ncol(B))
  ## now add remaining links to make number of edges match original network
  ER_B <- as.vector(ER_B)
  ER_B[sample(which(ER_B == 0), sum(B) - sum(ER_B), replace=FALSE)] <- 1
  ER_B <- matrix(ER_B, nrow(B), ncol(B))
  dimnames(ER_B) <- dimnames(B)
  ## convert back to igraph object
  return(graph_from_incidence_matrix(ER_B))
}

sparse_configuration_model <- function(B) {
    ## sort the matrix
    B <- B[order(rowSums(B)), order(colSums(B))]
    ## create a random spanning tree to ensure connectedness
    compatible <- FALSE
    num_tries_spanning <- 0
    max_tries_spanning <- 100
    while (!compatible) {
        num_tries_spanning <- num_tries_spanning + 1
        ## if already tried max_tries times, give up with warning
        if (num_tries_spanning > max_tries_spanning) {
            warning("Reached ", max_tries_spanning, " attempts to find a spanning tree without success")
            return(NULL)
        }
        ## find a random spanning tree
        spanning_tree <- random_spanning_tree(nrow(B), ncol(B))
        ## sort this matrix too
        spanning_tree <- spanning_tree[order(rowSums(spanning_tree)),
                                       order(colSums(spanning_tree))]
        ## check if the degree distributions are compatible with one another
        if (all(rowSums(spanning_tree) <= rowSums(B)) &
            all(colSums(spanning_tree) <= colSums(B))) {
            compatible <- TRUE
        }
    }
    ## now remove the necessary degrees from the matrix to be randomized
    restart <- TRUE
    num_tries_overall <- 0
    max_tries_overall <- 100
    while (restart) {
        restart <- FALSE
        num_tries_overall <- num_tries_overall + 1
        ## if already tried max_tries times, give up with Warning
        if (num_tries_overall > max_tries_overall) {
            warning("Reached ", max_tries_overall, " attempts to remove links without success")
            return(NULL)
        }
        B_minus_tree <- B
        row_degree_to_remove <- rowSums(spanning_tree)
        col_degree_to_remove <- colSums(spanning_tree)
        ## go through each row, removing the prescribed number of links from each
        for (rr in 1:nrow(B)) {
            ## how many links to remove in this row?
            links_to_remove <- row_degree_to_remove[rr]
            ## what are the options for removal?
            cols_to_choose_from <- which(col_degree_to_remove > 0 &
                                             B_minus_tree[rr,] == 1)
            ## pick an acceptable one
            violations <- TRUE
            num_tries_removal <- 0
            max_tries_removal <- 100
            while (violations) {
                num_tries_removal <- num_tries_removal + 1
                ## if already tried max_tries times, give up with Warning
                if (num_tries_removal > max_tries_removal) {
                    warning("Reached ", max_tries_removal, " attempts to remove links without success")
                    return(NULL)
                }
                ## if there are no acceptable links left to remove, give up and start over
                if (links_to_remove > length(cols_to_choose_from)) {
                    restart <- TRUE
                    cat("\n", num_tries_overall, "DEAD END", num_tries_removal, "/", max_tries_removal)
                    break()
                }
                ## otherwise, continue on
                to_remove <- my.sample(cols_to_choose_from, links_to_remove)
                testing_removal <- B_minus_tree
                testing_removal[rr, to_remove] <- 0
                ## I think this check is unnecessary because of how cols_to_choose_from is defined
                if (any(colSums(testing_removal) + colSums(spanning_tree) < colSums(B))) {
                    ## there is at least one violation
                    cat("\n", num_tries_overall, "VIOLATION", num_tries_removal, "/", max_tries_removal)
                    next()
                } else {
                    ## this is an acceptable removal
                    violations <- FALSE
                    B_minus_tree <- testing_removal
                    row_degree_to_remove <- rowSums(B_minus_tree) + rowSums(spanning_tree) - rowSums(B)
                    col_degree_to_remove <- colSums(B_minus_tree) + colSums(spanning_tree) - colSums(B)
                }
            }
            if (restart) break
        }
    }
    ## randomize the sub-network
    CM_B_minus_tree <- curve_ball(B_minus_tree, 25)
    ## add back in the spanning tree
    CM_B <- CM_B_minus_tree + spanning_tree
    ## check degree distribution
    if (!all(rowSums(B) == rowSums(CM_B))) {
        stop("degree distributions do not match")
    }
    ## convert back to igraph object
    return(graph_from_incidence_matrix(CM_B))
}
