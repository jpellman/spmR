#' @name spm_resss
#' @title SPM: Create residual sum of squares image (ResSS)
#' @usage Vo = spm_resss(Vi,Vo,R,flags)
#' @param Vi vector of mapped image volumes to work on (from spm_vol)
#' @param Vo handle structure for mapped output image volume
#' @param R residual forming matrix
#' flags       - 'm' for implicit zero masking
#' @return handle structure of output image volume after modifications
#'                 for writing
#

spm_resss <- function(Vi, R, flags="") {

    # -Argument checks
    if(flags == "m") {
        mask <- TRUE
    } else {
        mask <- FALSE
    }

    ni <- dim(R)[2]
    if(ni != prod(length(Vi))) {
        stop("incompatible dimensions (R) and (Vi)")
    }

    # Image dimension, orientation and voxel size checks
    # TODO

    # - C O M P U T A T I O N
    Vo <- array(NA, dim=dim(Vi[[1]]))  ## FIXME NA or 0?
    nB <- length(Vi)   
    # -Loop over planes computing ResSS
    for(p in 1:dim(Vo)[3]) {
        beta <- matrix(0, nrow=nB, ncol=prod(dim(Vo)[1:2]))
        for(b in 1:nB) {
            beta[b,] <- as.numeric(Vi[[b]][,,p])
        }
        ss <- colSums((R %*% beta)^2)
        Vo[,,p] <- ss 
    }

    Vo
}

