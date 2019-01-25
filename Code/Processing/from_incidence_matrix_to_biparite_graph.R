library(igraph)

from_incidence_matrix_to_bipartite <- function(tmp){
  g <- graph_from_incidence_matrix(tmp)
  cl <- clusters(g, "strong")
  g2 <- induced.subgraph(g, which(cl$membership == which.max(cl$csize)))
  return(g2)
}