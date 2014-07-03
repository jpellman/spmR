#' @name spm_dctmtx
#' @title SPM: Creates basis functions for Discrete Cosine Transform.
#' @usage C = spm_dctmtx(N,K,n)
#' @usage C = spm_dctmtx(N,K)
#' @usage D = spm_dctmtx(N,K,n,'diff')
#' @usage D = spm_dctmtx(N,K,'diff')
#' @param N dimension
#' @param K order
#' @param n optional points to sample
#____________________________________________________________________________
#' @description spm_dctmtx creates a matrix for the first few basis functions of a one
#' dimensional discrete cosine transform.
#' With the 'diff' argument, spm_dctmtx produces the derivatives of the
#' DCT.
#
#' @references Fundamentals of Digital Image Processing (p 150-154).  Anil K. Jain 1989.

spm_dctmtx <- function(N, K=N, n.opt=NULL, f.diff=NULL) {

    d <- 0
    n <- 0:(N-1)

    if(!is.null(n.opt)) {
        n <- n.opt
    }
    
    if(!is.null(f.diff)) {
        if(f.diff == "diff") {
            d = 1
        } else if(f.diff == "diff2") {
            d = 2
        } else {
            stop("unknown option for f.diff")
        }
    }

    C <- matrix(0, nrow=N, ncol=K)

    if(d == 0) {
        C[,1] <- 1/sqrt(N)
        for(k in 2:K) {
            C[,k] = sqrt(2/N)*cos(pi*(2*n+1)*(k-1)/(2*N))
        }
    } else if(d == 1) {
        ## FIXME: UNTESTED!!
        for(k in 2:K) {
            C[,k] = -2^(1/2)*(1/N)^(1/2)*sin(1/2*pi*(2*n*k-2*n+k-1)/N)*pi*(k-1)/N
        }
    } else if(d == 2) {
        ## FIXME: UNTESTED!!
        for(k in 2:K) {
            C[,k] = -2^(1/2)*(1/N)^(1/2)*cos(1/2*pi*(2*n+1)*(k-1)/N)*pi^2*(k-1)^2/N^2
        }
    }

    C
}

