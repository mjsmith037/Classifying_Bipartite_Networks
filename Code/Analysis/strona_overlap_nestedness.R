OverlapNestedness <- function(B) {
    RS <- rowSums(B)
    CS <- colSums(B)
    N.olap.vect <- vector(mode='numeric', length=sum(choose(dim(B), 2)))
    nn <- 1
    for (ii in 1:(dim(B)[1] - 1)) {
        for (jj in (ii + 1):dim(B)[1]) {
            S <- sum(B[ii,] + B[jj,] == 2)
            k <- 1:min(RS[ii], RS[jj])
            P <- sum((choose(dim(B)[2], k) *
                          choose(dim(B)[2] - k, RS[jj] - k) *
                          choose(dim(B)[2] - RS[jj], RS[ii] - k)) /
                         (choose(dim(B)[2], RS[jj]) *
                              choose(dim(B)[2], RS[ii])) * k)
            omega <- (S - P) / min(RS[ii], RS[jj])
            if (is.nan(P)) return(NaN)
            if (S < P) {
                if (RS[ii] + RS[jj] - dim(B)[2] < 0) {
                    OMEGA <- (min(RS[ii], RS[jj]) - P) / min(RS[ii], RS[jj])
                } else OMEGA <- (P - (RS[ii] + RS[jj] - dim(B)[2])) / min(RS[ii], RS[jj])
            } else if (S > P) {
                OMEGA <- (min(RS[ii], RS[jj]) - P) / min(RS[ii], RS[jj])
            } else OMEGA <- 1
            N.olap.vect[nn] <- omega / OMEGA
            nn <- nn + 1
        }
    }
    for (ii in 1:(dim(B)[2] - 1)) {
        for (jj in (ii + 1):dim(B)[2]) {
            S <- sum(B[,ii] + B[,jj] == 2)
            k <- 1:min(CS[ii], CS[jj])
            P <- sum((choose(dim(B)[1], k) *
                          choose(dim(B)[1] - k, CS[jj] - k) *
                          choose(dim(B)[1] - CS[jj], CS[ii] - k)) /
                         (choose(dim(B)[1], CS[jj]) *
                              choose(dim(B)[1], CS[ii])) * k)
            omega <- (S - P) / min(CS[ii], CS[jj])
            if (is.nan(P)) return(NaN)
            if (S < P) {
                if (CS[ii] + CS[jj] - dim(B)[1] < 0) {
                    OMEGA <- (min(CS[ii], CS[jj]) - P) / min(CS[ii], CS[jj])
                } else OMEGA <- (P - (CS[ii] + CS[jj] - dim(B)[1])) / min(CS[ii], CS[jj])
            } else if (S > P) {
                OMEGA <- (min(CS[ii], CS[jj]) - P) / min(CS[ii], CS[jj])
            } else OMEGA <- 1
            N.olap.vect[nn] <- omega / OMEGA
            nn <- nn + 1
        }
    }
    return(mean(N.olap.vect))
}
