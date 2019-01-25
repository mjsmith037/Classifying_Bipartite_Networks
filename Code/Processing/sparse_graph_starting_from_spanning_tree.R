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