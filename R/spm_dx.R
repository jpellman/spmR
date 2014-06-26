#' @name spm_dx
#' @title SPM: returns dx(t) = (expm(dfdx*t) - I)*inv(dfdx)*f
#' @usage [dx] = spm_dx(dfdx,f,[t])
#' @param dfdx df/dx
#' @param f dx/dt
#' @param t integration time: (default t = Inf);
#'          if t is a cell (i.e., {t}) then t is set to:
#'          exp(t - log(diag(-dfdx))
#
#' @return x(t) - x(0)
#--------------------------------------------------------------------------
#' @description Integration of a dynamic system using local linearization.  This scheme
#' accommodates nonlinearities in the state equation by using a functional of
#' f(x) = dx/dt.  This uses the equality
#'
#'             expm([0   0     ]) = (expm(t*dfdx) - I)*inv(dfdx)*f
#'                  [t*f t*dfdx]
#'
#' When t -> Inf this reduces to
#'
#'              dx(t) = -inv(dfdx)*f
#'
#' These are the solutions to the gradient ascent ODE
#'
#'            dx/dt   = k*f = k*dfdx*x =>
#'
#'            dx(t)   = expm(t*k*dfdx)*x(0)
#'                    = expm(t*k*dfdx)*inv(dfdx)*f(0) -
#'                      expm(0*k*dfdx)*inv(dfdx)*f(0)
#'
#' When f = dF/dx (and dfdx = dF/dxdx), dx represents the update from a
#' Gauss-Newton ascent on F.  This can be regularised by specifying {t}
#' A heavy regularization corresponds to t = -4 and a light
#' regularization would be t = 4. This version of spm_dx uses an augmented
#' system and the Pade approximation to compute requisite matrix
#' exponentials

spm_dx <- function(dfdx, f, t=Inf) {

    # t is a regulariser
    if( is.list(t) ) {
        t <- exp( unlist(t) - log( diag(-dfdx)) )
    }

    cat("## dfdx ="); print(dfdx); cat("\n")
    cat("## f = "); print(f); cat("\n")
    cat("## t = "); print(t); cat("\n")

    # use a [pseudo]inverse if all t > TOL
    if(min(t) > exp(16)) {
        ## TODO!!
        dx = -spm_pinv(dfdx)*spm_vec(f);
        dx =  spm_unvec(dx,f);
    } else {
        # ensure t is a scalar or matrix
        if(length(t) > 1) {
            tval <- t
            t <- diag( length(tval) )
            diag(t) <- tval
        }

        q <- matrix(0, nrow=max(dim(dfdx))+1, ncol=1); q[1,1] <- 1

        # augment Jacobian and take matrix exponential
        Jx <- rbind(0,cbind(t %*% f, t %*% dfdx))
        dx <- Matrix::expm(Jx) %*% q
        dx <- dx[2:nrow(dx),]  ### why is there 'f' in spm_dx.m??
    }

    dx

}

