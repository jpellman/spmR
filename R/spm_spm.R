#' @name spm_spm
#' @title SPM: [Re]ML Estimation of a General Linear Model
#' @usage FORMAT [SPM] = spm_spm(SPM)
#
#' @description Required fields of SPM:
#'
#' xY.VY - nScan x 1 struct array of image handles (see spm_vol)
#'         Images must have the same orientation, voxel size and data type
#'       - Any scaling should have already been applied via the image handle
#'         scalefactors.
#'
#' xX    - Structure containing design matrix information
#'       - Required fields are:
#'         xX.X      - Design matrix (raw, not temporally smoothed)
#'         xX.name   - cellstr of parameter names corresponding to columns
#'                     of design matrix
#'       - Optional fields are:
#'         xX.K      - cell of session-specific structures (see spm_filter)
#'                   - Design & data are pre-multiplied by K
#'                     (K*Y = K*X*beta + K*e)
#'                   - Note that K should not smooth across block boundaries
#'                   - defaults to speye(size(xX.X,1))
#'         xX.W      - Optional whitening/weighting matrix used to give
#'                     weighted least squares estimates (WLS). If not specified
#'                     spm_spm will set this to whiten the data and render
#'                     the OLS estimates maximum likelihood
#'                     i.e. W*W' = inv(xVi.V).
#'
#' xVi   - Structure describing intrinsic temporal non-sphericity
#'       - Required fields are:
#'         xVi.Vi    - array of non-sphericity components
#'                   - defaults to {speye(size(xX.X,1))} - i.i.d.
#'                   - specifying a cell array of constraints (Qi)
#'                     These constraints invoke spm_reml to estimate
#'                     hyperparameters assuming V is constant over voxels.
#'                     that provide a high precise estimate of xX.V
#'       - Optional fields are:
#'         xX.V      - Optional non-sphericity matrix.  Cov(e) = sigma^2*V
#'                     If not specified spm_spm will compute this using
#'                     a 1st pass to identify significant voxels over which
#'                     to estimate V.  A 2nd pass is then used to re-estimate
#'                     the parameters with WLS and save the ML estimates
#'                     (unless xX.W is already specified).
#'
#' xM    - Structure containing masking information, or a simple column vector
#'         of thresholds corresponding to the images in VY [default: -Inf]
#'       - If a structure, the required fields are:
#'         xM.TH - nVar x nScan matrix of analysis thresholds, one per image
#'         xM.I  - Implicit masking (0=>none, 1 => implicit zero/NaN mask)
#'         xM.VM - struct array of explicit mask image handles
#'       - (empty if no explicit masks)
#'               - Explicit mask images are >0 for valid voxels to assess.
#'               - Mask images can have any orientation, voxel size or data
#'                 type. They are interpolated using nearest neighbour
#'                 interpolation to the voxel locations of the data Y.
#'       - Note that voxels with constant data (i.e. the same value across
#'         scans) are also automatically masked out.
#'
#' swd   - Directory where the output files will be saved [default: pwd]
#'         If exists, it becomes the current working directory.
#'
#' In addition, global SPM "defaults" variable is used (see spm_defaults):
#' 
#' stats.<modality>.UFp - critical F-threshold for selecting voxels over 
#'                        which the non-sphericity is estimated (if 
#'                        required) [default: 0.001]
#' 
#' stats.maxres         - maximum number of residual images for smoothness
#'                        estimation
#'
#' stats.maxmem         - maximum amount of data processed at a time (in bytes)
#'
#' modality             - SPM modality {'PET','FMRI','EEG'}
#
#' @param SPM
#' @return SPM

spm_spm <- function(SPM) {

#==========================================================================
# - A N A L Y S I S   P R E L I M I N A R I E S
#%==========================================================================

    # initialize
    xX <- SPM$xX
    nScan <- nrow(xX$X)
    nBeta <- ncol(xX$X)

    # xM
    xM <- SPM$xM

    # check confounds (xX$K) and non-sphericity (xVi)
    if(!is.list(xX$K)) {
       xX$K <- 1
    }

    xVi <- SPM$xVi


   
    # Get non-sphericity V
    if(!is.null(SPM$xVi$V)) {
        V <- SPM$xVi$V
    } else {
        V <- NULL
    }

    # compute weight matrix 'W' [ WW' = inv(V) ]
    if(!is.null(V)) {
        done <- TRUE

        # compute W
        # following spm, we use the SVD route (for now)
        out <- svd(V)
        s <- diag( 1/sqrt(out$d) )
        W <- out$u %*% s %*% t(out$u)
        idx <- which(abs(W) < 1e-6)
        W[idx] <- 0.0
        xX$W <- W

    } else {
        # require 2nd pass
        done <- FALSE
        W <- diag(nScan)
    }

    # design space and projector matrix (pseudoinverse) for WLS
    xX$xKXs <- spm_sp("Set", spm_filter(xX$K, W %*% xX$X)) 
    xX$pKX  <- spm_sp("x-", xX$xKXs)
    erdf    <- spm_SpUtil("trRV", xX$xKXs)


    # if xVi$V is not defined, compute Hsqr and F-threshold under iid 
    if(!done) {
        Fcname <- "effects of interest"
        iX0  <- c(SPM$xX$iB, SPM$xX$iG) # effects that need *not* to be tested
        xCon <- spm_FcUtil("set", name=Fcname, STAT="F",
                           set_action="iX0", value=iX0, sX=xX$xKXs)
        X1o  <- spm_FcUtil("X1o",  Fc=xCon, sX=xX$xKXs)
        Hsqr <- spm_FcUtil("hsqr", Fc=xCon, sX=xX$xKXs)
        trRV <- spm_SpUtil("trRV", xX$xKXs)
        trMV <- spm_SpUtil("trMV", X1o)

        UFp  <- defaults.stats.fmri.ufp

        UF  <- qf(UFp, trMV, trRV, lower.tail=FALSE)
    }


    # Image dimensions and data
    DIM <- dim(SPM$xY$VY)[1:3]
    xdim <- DIM[1]; ydim <- DIM[2]; zdim <- DIM[3]

    # Maximum number of residual images for smoothness estimation
    MAXRES <- defaults.stats.maxres
    nSres  <- min(nScan, MAXRES)


    # Initialise output images (unless this is a 1st pass for ReML)
    if(done) {
        # mask (integer array)
        VM <- array(0L, dim=DIM)

        # beta's
        Vbeta <- vector("list", length=nBeta)
        for(b in 1:nBeta) {
            Vbeta[[b]] <- array(NA, dim=DIM)
        }

        # residual SS
        VResMS <- array(NA, dim=DIM)

        # standardised residual images
        VresI <- vector("list", length=nSres)
        for(i in 1:nSres) {
            VresI[[i]] <- array(NA, dim=DIM)
        }
    }

#==========================================================================
# - F I T   M O D E L   &   W R I T E   P A R A M E T E R    I M A G E S
#==========================================================================

    # Note: in spmR, we do not split the analysis in blocks
    # we estimate one plane at a time

    # do GLM per plane
    CY <- matrix(0, nrow=nScan, ncol=nScan) # (Y - <Y>) * (Y - <Y>)'
    Cy <- matrix(0, nrow=nScan, ncol=nScan) # Y*Y' spatially whitened
    EY <- numeric(nScan)                    # Y    for ReML
    Q <- 0                                  # number of inmask voxels
    s <- 0                                  # Volume (voxels > UF)
    i_res <- floor(seq(1,nScan,length=nSres) + 0.5) # Indices for residual

    if(SPM$SPMR$verbose) cat("Estimating parameters:\n")
    for(z in 1:zdim) {

        if(SPM$SPMR$verbose) cat("Plane: ", z, sep="")

        # select plane+ts, make it flat (voxel + ts)
        YY <- matrix(SPM$xY$VY[,,z,], nrow=nScan, byrow=TRUE)

        # select 'in mask' voxels, above threshold, non-constant voxels
        Cm <- rep(TRUE, ncol(YY))  # all voxels in XY plane
            # FIXME:
            # 1. explicit mask? handle here
            # 

            # 2. > threshold (per scan)
            for(i in 1:nScan) {
                idx <- which(YY[i,] < SPM$xM$TH[i])
                Cm[idx] <- FALSE
            }
            # 3. implicit mask ??
            # 4. constant voxels
            idx <- which(apply(YY[,Cm,drop=FALSE], MARGIN=2,
                              function(x) { all(x[1] == x) }))
            if(length(idx) > 0 && SPM$SPMR$verbose) {
                cat(" Constant voxels removed: ", length(idx), sep="")
            }
            Cm[Cm][idx] <- FALSE

        YY <- YY[,Cm,drop=FALSE]
        S <- sum(Cm)
        Q <- Q + S
        if(done) {
            VM[,,z][Cm] <- 1L
        }

        if(SPM$SPMR$verbose) cat(" Voxels selected: ", S, sep="")

        # proceed with General Linear Model
        if(S > 0) {
            # whiten/weight data and remove filter confounds
            KWY <- spm_filter(xX$K, W %*% YY)

            # General linear model: weighted least squares estimation
            beta  <- xX$pKX %*% KWY           # beta = pinv(X) %*% Y
            res   <- spm_sp("r",xX$xKXs,KWY)  # res  = Y - Y.hat
            ResSS <- colSums(res^2)           # residual SS
            rm(KWY)

            # If ReML hyperparameters are needed for xVi.V
            if(!done) {
                # first pass: select 'active' voxels
                F <- (colSums((Hsqr %*% beta)^2) / trMV) / (ResSS/trRV)
                idx <- which(F > UF)
                q <- length(idx)
                cat(" (> UF: ", q, ")", sep="")
                if(q > 0) {
                   s <- s + q
                   tmp <- sqrt(trRV/ResSS[idx])
                   q <- diag(tmp, nrow=length(tmp))
                   YY <- YY[,idx] %*% q
                   #Cy <- Cy + YY %*% t(YY)
                   ### FIXME:
                   ### the end result (Cy) is slightly off!!
                   ### accumulating round-off error??
                   Cy <- Cy + tcrossprod(YY)
                }
            }

            # if we are saving the WLS parameters
            if(done) {
                # beta
                for(b in 1:nBeta) {
                    Vbeta[[b]][,,z][Cm] <- beta[b,]
                }
                # SSE
                VResMS[,,z][Cm] <- ResSS
                # standardized residuals
                for(i in 1:nSres) {
                    VresI[[i]][,,z][Cm] <- 
                        (res[i_res[i],,drop=FALSE] / sqrt(ResSS/erdf))
                }
            } 

        } # S

        # save parameters
        if(done) {
            CY <- CY + YY %*% t(YY)
            EY <- EY + rowSums(YY)
        }

        if(SPM$SPMR$verbose) cat("...done.\n")

    } # z-plane
    if(SPM$SPMR$verbose) cat("done\n")

    # average sample covariance and mean of Y (over voxels)
    if(done) {
        CY <- CY/Q
        EY <- EY/Q
        CY <- CY - EY %*% t(EY)
    } else {
        stopifnot(s > 0)
        # do REML and start over
        if(SPM$SPMR$verbose) cat("Temporal non-sphericity (over voxels) ... ReML estimation\n")
        Cy <- Cy/s
        if(is.list(SPM$xX$K)) {
            ## !!
            ## FIXME: we assume only one Hparam/K per session!
            ## !!
            Xp <- cbind(xX$X, xX$K[[1]]$X0)
            V <- spm_reml(Cy, Xp, SPM$xVi$Vi, verbose=SPM$SPMR$verbose)
            h <- attr(V, "h")
        } else {
            V <- spm_reml(Cy, xX$X, xVi$Vi, verbose=SPM$SPMR$verbose)
            h <- attr(V, "h")
        }

        # normalize non-sphericity and save hyperparameters
        V <- V*nScan/sum(diag(V))
        xVi$h  <- h
        xVi$V  <- V
        xVi$Cy <- Cy
        SPM$xVi <- xVi

        # and again
        out <- spm_spm(SPM)
        return(out)
    }

    # Use non-sphericity xVi.V to compute [effective] degrees of freedom
    xX$V <- spm_filter(xX$K, spm_filter(xX$K, W %*% V %*% t(W)))  # KWVW'K'
    # WVW   <- W %*% SPM$xVi$V %*% t(W)
    # KWVW  <- t( WVW - K %*% (t(K) %*% WVW) )
    # KWVWK <- KWVW - K %*% (t(K) %*% KWVW)
    trRV <- spm_SpUtil("trRV", xX$xKXs, xX$V)
    trRVRV <- attr(trRV, "trRVRV")
    xX$trRV <- trRV
    xX$trRVRV <- trRVRV
    xX$erdf   <- trRV^2/trRVRV
    xX$Bcov <- xX$pKX %*% xX$V %*% t(xX$pKX)

    # Set VResMS scalefactor as 1/trRV (raw voxel data is ResSS)
    VResMS <- VResMS * 1/xX$trRV

    # Smoothness estimates of component fields and RESEL counts for volume
    if(SPM$SPMR$verbose)
        cat("\nEstimating smoothness + RESEL counts ...")
    out <- spm_est_smoothness(VresI, VM, c(nScan, erdf))
    FWHM <- out$FWHM
    VRpv <- out$VRpv
    R    <- out$R
    if(SPM$SPMR$verbose) {
        cat("done.\n")
        cat("FWHM = "); print(FWHM);
        cat("RESEL count R = "); print(R)
        cat("\n")
    }

    # Delete the residuals images
    # (just don't write them to disk)

    # Compute scaled design matrix for display purposes
    ## xX$nKX        = spm_DesMtx("sca", xX$xKXs$X, xX$name)

    # place fields in SPM
    XYZ <- which(VM == 1L, arr.ind=TRUE)
    XYZ <- unname(XYZ)

    SPM$xVol <- list()
    SPM$xVol$XYZ   <- XYZ                 #-InMask XYZ coords (voxels)
    #SPM$xVol$M     <- M                  #-voxels -> mm
    #SPM$xVol$iM    <- inv(M);             #-mm -> voxels
    SPM$xVol$DIM   <- DIM                 #-image dimensions
    SPM$xVol$FWHM  <- FWHM;               #-Smoothness data
    SPM$xVol$R     <- R;                  #-Resel counts
    SPM$xVol$S     <- Q;                  #-Volume (voxels)
    SPM$xVol$VRpv  <- VRpv;               #-Filehandle - Resels per voxel

    SPM$Vbeta      <- Vbeta               #-Filehandle - Beta
    SPM$VResMS     <- VResMS              #-Filehandle - Hyperparameter
    SPM$VM         <- VM                  #-Filehandle - Mask
    SPM$VresI      <- VresI
  
    SPM$xVi        <- xVi                 # non-sphericity structure
    SPM$xVi.CY     <- CY                  #-<(Y - <Y>)*(Y - <Y>)'>

    SPM$xX         <- xX                  #-design structure
 
    SPM$xM         <- xM                  #-mask structure

    SPM$xCon       <- list()              #-contrast structure
 
    #SPM$SPMid      = SPMid;
    #SPM$swd        = pwd;


    SPM
}

