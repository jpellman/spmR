#' @name spm_reml
#' @title SPM: ReML estimation of [improper] covariance components from y*y'
#' @usage [C,h,Ph,F,Fa,Fc] = spm_reml(YY,X,Q,N,D,t);
#
#' @param YY (m x m) sample covariance matrix Y*Y'  {Y = (m x N) data matrix}
#' @param X (m x p) design matrix
#' @param Q {1 x q} covariance components
#' @param N number of samples
#' @param D Flag for positive-definite scheme
#' @param t regularisation (default 4)
#
#' @return A list containing:
#' C   - (m x m) estimated errors = h(1)*Q{1} + h(2)*Q{2} + ..., 
#' h   - (q x 1) ReML hyperparameters h, 
#' Ph  - (q x q) conditional precision of h, 
#' F   - [-ve] free energy F = log evidence = p(Y|X,Q) = ReML objective, 
#' Fa  - accuracy, 
#' Fc  - complexity (F = Fa - Fc)
#
#' @description Performs a Fisher-Scoring ascent on F to find ReML variance parameter
#' estimates.

spm_reml <- function(YY, X, Q, N=1, D=0, t=4, verbose=FALSE) {

    # check defaults
    K = 32 # number of iterations

    # catch NaNs
    W <- Q
    idx <- which(is.na(YY))
    if(length(idx) > 0) {
        stop("NaNs found in YY; see spm_reml.m code to see what to do")
    }

    # dimensions
    n <- ncol(Q[[1]])
    m <- length(Q)

    # ortho-normalise X
    X <- svd(X)$u


    # initialize and specify hyperpriors
    h <- numeric(m)
    for(i in 1:m) {
        h[i] = sum(any(diag(Q[[i]]) != 0))
    }
    hE <- numeric(m)
    hP <- diag(m)/exp(32)
    dF <- Inf
    D  <- 8*(D > 0)

    # spmr added
    dh <- numeric(m)
    dFdh <- numeric(m)
    dFdhh <- matrix(0, m, m)


    # ReML (EM/VB)
    for(k in 1:K) {

        # compute current estimate of covariance
        C <- matrix(0, n, n)
        for(i in 1:m) {
            C <- C + Q[[i]] * h[i]
        }

        # positive [semi]-definite check
        # TODO!!

        # E-step: conditional covariance cov(B|y) {Cq}
        iC <- solve(C + diag(n)/exp(32))
        iCX <- iC %*% X
        Cq  <- solve(t(X) %*% iCX)

        # M-step: ReML estimate of hyperparameters

        # Gradient dF/dh (first derivatives)
        P <- iC - iCX %*% Cq %*% t(iCX)
        U <- diag(n) - (P %*% YY)/N
        PQ <- vector("list", length=m)
        for(i in 1:m) {
            # dF/dh = -trace(dF/diC*iC*Q{i}*iC)
            PQ[[i]] <- P %*% Q[[i]]
            dFdh[i] <- -1 * sum(diag(PQ[[i]] %*% U)) * N/2
        }

        # Expected curvature E{dF/dhh} (second derivatives)
        for(i in 1:m) {
            for(j in i:m) {
                # dF/dhh = -trace{P*Q{i}*P*Q{j}}
                dFdhh[i,j] <- -1 * sum(diag( PQ[[i]] %*% PQ[[j]])) * N/2
                dFdhh[j,i] <- dFdhh[i,j]
            }
        }

        # add hyperpriors
        e     <- h     - hE
        dFdh  <- dFdh  - hP %*% e
        dFdhh <- dFdhh - hP

        # Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
        dh <- spm_dx(dFdhh, dFdh, as.list(t))
        h  <- h + dh

        # predicted change in F - increase regularisation if increasing
        pF <- t(dFdh) %*% dh
        if(pF > dF) {
            t <- t - 1
        } else {
            t <- t + 1/4
        } 
        dF <- pF

        if(verbose) {
            cat("  ReML Iteration ", k, "pF = ", pF, "\n")
        }

        # check for convergence
        if(dF < 1e-1) break

    } # K REML iterations

    # re-build predicted covariance
    V <- matrix(0, n, n)
    for(i in 1:m) {
        V <- V + W[[i]] * h[i]
    }

    # check V is positive semi-definite (if not already checked)
    ## TODO!!

    # addition out arguments
    attr(V, "h") <- h

    ## TODO!!

    V
}



