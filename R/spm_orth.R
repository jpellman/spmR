# recursive Gram-Schmidt orthogonalisation of basis functions
# FORMAT X = spm_orth(X,OPT)
#
# X   - matrix
# OPT - 'norm' for Euclidean normalisation
#     - 'pad'  for zero padding of null space [default]
#
# serial orthogonalisation starting with the first column


## TODO: OPT="pad" not used; no checking!

spm_orth <- function(X, OPT="pad") {

    # should be done using QR...
    X <- as.matrix(X)
    nc <- ncol(X)

    if(nc == 1) {
        return(X)
    }

    x <- X[,1]
    for(i in 2:nc) {
        x <- X[,1:(i-1)]
        X[,i] <- X[,i] - x %*% (solve(t(x) %*% x) %*% t(x) %*% X[,i])
    }

    X

}
