#' @name spm_est_smoothness
#' @title SPM: Estimation of smoothness based on [residual] images
#' @usage [FWHM,VRpv,R] = spm_est_smoothness(V,VM,[ndf]);
#
#' @param V Filenames or mapped standardized residual images
#' @param VM Filename of mapped mask image
#' @param ndf A 2-vector, [n df], the original n & dof of the linear model
#
#' @return A list containing: FWHM  - estimated FWHM in all image directions,
#' VRpv  - handle of Resels per Voxel image,
#' R     - vector of resel counts

# Matlab version: Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

spm_est_smoothness <- function(V, VM, ndf) {

    #-Initialise
    ##--------------------------------------------------------------------------
    if(any(is.na(ndf))) {
        ndf <- c(length(V), length(V))
    }
    n_full <- ndf[1]
    edf    <- ndf[2]

    #-Initialise RESELS per voxel image
    ##--------------------------------------------------------------------------
    DIM <- dim(V[[1]])
    VRpv <- array(NA, dim=DIM)

    #-Dimensionality of image
    ##--------------------------------------------------------------------------
    D <- length(DIM)
    if(D == 0L) {
        return(list(FWHM=c(Inf,Inf,Inf), VRpv=VRpv, R=c(0,0,0)))
    }

    #-Find voxels within mask
    ##--------------------------------------------------------------------------
    d = VM
    Ixyz <- which(VM == 1, arr.ind=TRUE)
    Ix <- Ixyz[,1]; Iy <- Ixyz[,2]; Iz <- Ixyz[,3]

    #-Compute covariance of derivatives
    ##--------------------------------------------------------------------------
    #nonzero_in_mask <- length(which(VM > 0L))
    nonzero_in_mask <- length(Ix)
    L <- array(0, dim=c(nonzero_in_mask, D, D))
    ssq <- numeric(nonzero_in_mask)
    for (i in 1:length(V)) {
        #cat("\nspm_sample_vol (computing derivatives) for V = ", i, "\n")
        #[d,dx,dy,dz] = spm_sample_vol(V(i),Ix,Iy,Iz,hold=1,gradient=TRUE);
        out <- spm_sample_vol(V[[i]],Ix,Iy,Iz,hold=1,gradient=TRUE)
        d  <- out$out[as.logical(VM)]
        dx <- out$gradx
        dy <- out$grady
        dz <- out$gradz
    
        # sum of squares
        #----------------------------------------------------------------------
        ssq  = ssq + d^2
    
        # covariance of finite differences
        #----------------------------------------------------------------------
        if (D >= 1L) {
            L[,1,1] = L[,1,1] + dx*dx
        }
        if (D >= 2L) {
            L[,1,2] = L[,1,2] + dx*dy
            L[,2,2] = L[,2,2] + dy*dy
        }
        if (D >= 3) {
            L[,1,3] = L[,1,3] + dx*dz
            L[,2,3] = L[,2,3] + dy*dz
            L[,3,3] = L[,3,3] + dz*dz
        }
    }   

    L  = L/length(V);     # Average
    L  = L*(n_full/edf);  # Scale


    #-Evaluate determinant (and xyz components for FWHM)
    #--------------------------------------------------------------------------
    if (D == 1) {
        resel_xyz <- L
        resel_img <- L
    }
    if (D == 2) {
        resel_xyz <- cbind(L[,1,1], L[,2,2])
        resel_img <- (L[,1,1]*L[,2,2] - L[,1,2]*L[,1,2])
    }
    if (D == 3) {
        resel_xyz <- cbind(L[,1,1], L[,2,2], L[,3,3])
        resel_img <- (L[,1,1]*L[,2,2]*L[,3,3] + 
                      L[,1,2]*L[,2,3]*L[,1,3]*2 - 
                      L[,1,1]*L[,2,3]*L[,2,3] - 
                      L[,1,2]*L[,1,2]*L[,3,3] -
                      L[,1,3]*L[,2,2]*L[,1,3])
    }
    resel_img[resel_img<0] <- 0
    # Convert det(Lambda) and diag(Lambda) to units of resels
    resel_img <- sqrt(resel_img/(4*log(2))^D)
    resel_xyz <- sqrt(resel_xyz/(4*log(2)))

    #-Write Resels Per Voxel image
    #--------------------------------------------------------------------------
    VRpv[as.logical(VM)] <- resel_img

    #-(unbiased) RESEL estimator and Global equivalent FWHM
    # where we desire FWHM with components proportional to 1./mean(resel_xyz),
    # but scaled so prod(1./FWHM) agrees with (the unbiased) mean(resel_img).
    #--------------------------------------------------------------------------
    i     <- which(is.na(ssq) | ssq < sqrt(.Machine$double.eps))
    resel_img <- mean(resel_img[-i])
    resel_xyz <- colMeans(resel_xyz[-i,])

    RESEL <- resel_img^(1/D)*(resel_xyz/prod(resel_xyz)^(1/D))
    #FWHM  = full(sparse(1,1:D,1./RESEL,1,3));
    FWHM <- 1/RESEL
    FWHM[FWHM == 0.0] <- 1

    #-resel counts for search volume (defined by mask)
    #--------------------------------------------------------------------------
    # R0   = spm_resels_vol(VM,[1 1 1])';
    # R    = R0.*(resel.^([0:3]/3));
    # OR
    R      <- spm_resels_vol(VM,FWHM)

    return(list(FWHM=FWHM, VRpv=VRpv, R=R))
} 
