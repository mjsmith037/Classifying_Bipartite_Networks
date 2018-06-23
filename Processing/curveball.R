## curve ball algorithm code slightly modified from:
## Strona, G., Nappo, D., Boccacci, F., Fattorini, S., & San-Miguel-Ayanz, J.
##     (2014). A fast and unbiased procedure to randomize ecological binary
##     matrices with fixed row and column totals. Nature Communications, 5.
curve_ball <- function(mat, repmultiplier){
    ## convert the matrix to a list of link positions (by row)
    link_list <- lapply(1:nrow(mat), function(row) which(mat[row,] == 1))
    ## for each iteration of the scrambling routine:
    for (rr in 1:(repmultiplier * nrow(mat))){
        ## select two distinct rows at random
        rows <- sample(1:nrow(mat), 2)
        ## and get the positions of the links in those rows
        row_1_links <- link_list[[rows[1]]]
        row_2_links <- link_list[[rows[2]]]
        ## find which links are in common between the two rows
        links_in_common <- intersect(row_1_links, row_2_links)
        ## if the rows are not proper subsets of each other
        if (!(length(links_in_common) %in% c(length(row_1_links), length(row_2_links)))) {
            ## scramble the links that differ between the two rows
            links_different <- sample(setdiff(c(row_1_links, row_2_links), links_in_common))
            ## save how many links in the first row are not shared with the second row
            len_unique <- length(row_1_links) - length(links_in_common)
            ## update the list of link positions
            link_list[[rows[1]]] <- c(links_in_common, links_different[1:len_unique])
            link_list[[rows[2]]] <- c(links_in_common,
                               links_different[(len_unique + 1):length(links_different)])
        }
    }
    ## initialize a matrix for output (same size as input)
    rm <- matrix(0, nrow(mat), ncol(mat))
    ## convert the (now randomized) list of links back to a matrix
    lapply(1:nrow(mat), function(row) rm[row, link_list[[row]]] <<- 1)
    return(rm)
}
