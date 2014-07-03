#' @name spm_Q
#' @title SPM: returns an (n x n) autocorrelation matrix for an AR(p) process
#' @usage [Q] = spm_Q(a,n,q)
#
#' @param a a vector pf p AR coefficients
#' @param n size of Q
#' @param q switch to return precision [default q = 0]
#' @return an (n x n) autocorrelation matrix for an AR(p) process
#__________________________________________________________________________
#' @description spm_Q uses a Yule-Walker device to compute K where:
#' 
#' y = K*z
#' 
#' such that y is an AR(n) process generated from an i.i.d innovation 
#' z.  This means
#' 
#' cov(y) = <K*z*z'*K> = K*K'
#' 
#' Critically, this is not the correlation because if cov(z) = eye(n) 
#' then trace(cov(y)) ~= n.
#' The reason the diagonals of corr(y) are not constant is that we 
#' are modeling finite length AR sequences, which incur boundary effects 
#' at the beginning and end of the sequence.  These are resolved with a
#' Toeplitz form
#' @seealso spm_Ce

spm_Q <- function(a, n, q=0) {

    p <- length(a)
    if(p > 1) {
        stop("length(a) > 1 not implemented!")
    }

    if(q) {
        stop("q=0 is not implemented")
    } else {
        # compute Q
        P <- diag(n)
        idx <- which(P > 0)[-n]  # trick to set off-diagonal element
        P[idx+1] <- -a           # may not be very robust?
        K <- solve(P)
        K[ which(K < 1e-4) ] <- 0.0
        Q <- K %*% t(K)
        Q <- toeplitz(Q[,1])
    }

    Q
}

