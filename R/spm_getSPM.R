# Compute a specified and thresholded SPM/PPM following estimation
# FORMAT [SPM,xSPM] = spm_getSPM;
# Query SPM in interactive mode.
#
# FORMAT [SPM,xSPM] = spm_getSPM(xSPM);
# Query SPM in batch mode. See below for a description of fields that may
# be present in xSPM input. Values for missing fields will be queried
# interactively.
#
# xSPM      - structure containing SPM, distribution & filtering details
# .swd      - SPM working directory - directory containing current SPM.mat
# .title    - title for comparison (string)
# .Z        - minimum of Statistics {filtered on u and k}
# .n        - conjunction number <= number of contrasts
# .STAT     - distribution {Z, T, X, F or P}
# .df       - degrees of freedom [df{interest}, df{residual}]
# .STATstr  - description string
# .Ic       - indices of contrasts (in SPM.xCon)
# .Im       - indices of masking contrasts (in xCon)
# .pm       - p-value for masking (uncorrected)
# .Ex       - flag for exclusive or inclusive masking
# .u        - height threshold
# .k        - extent threshold {voxels}
# .XYZ      - location of voxels {voxel coords}
# .XYZmm    - location of voxels {mm}
# .S        - search Volume {voxels}
# .R        - search Volume {resels}
# .FWHM     - smoothness {voxels}
# .M        - voxels -> mm matrix
# .iM       - mm -> voxels matrix
# .VOX      - voxel dimensions {mm} - column vector
# .DIM      - image dimensions {voxels} - column vector
# .Vspm     - Mapped statistic image(s)
# .Ps       - uncorrected P values in searched volume (for voxel FDR)
# .Pp       - uncorrected P values of peaks (for peak FDR)
# .Pc       - uncorrected P values of cluster extents (for cluster FDR)
# .uc       - 0.05 critical thresholds for FWEp, FDRp, FWEc, FDRc
# .thresDesc - description of height threshold (string)
#
# Copyright MATLAB version: (C) 2008 Wellcome Trust Centre for Neuroimaging
 
spm_getSPM <- function(SPM, Ic=1L, thresDesc="FWE", u=0.05, k=0) {
    # Ic          contrast index
    # thresDesc   "FWE", "FDR", "none"
    # u           p-value for FWE
    # k           extent threshold (voxels)


    xX   <- SPM$xX                       #-Design definition structure
    XYZ  <- SPM$xVol$XYZ                 #-XYZ coordinates
    S    <- SPM$xVol$S                   #-search Volume {voxels}
    R    <- SPM$xVol$R                   #-search Volume {resels}
    M    <- SPM$xVol$M[1:3,1:3]          #-voxels to mm matrix
    VOX  <- sqrt(diag(tcrossprod(M)))    #-voxel dimensions

    #==========================================================================
    # - C O N T R A S T S ,   S P M    C O M P U T A T I O N ,    M A S K I N G
    #==========================================================================
    xCon <- SPM$xCon
  
    # FIXME: we allow no conjunctions (for now!) only a single contrast
    nc <- 1L; n <- 1L

    # apply MASKING --- TODO!!!
    Mask <- 0L
    pm <- numeric(0);
    Im <- numeric(0);
    Ex <- numeric(0);

    #-Create/Get title string for comparison
    #--------------------------------------------------------------------------
    if(nc == 1L) {
        str <- xCon[[Ic[1L]]]$name
    }
    titlestr <- str

    #-Bayesian or classical Inference?
    #==========================================================================
    
    #-Compute & store contrast parameters, contrast/ESS images, & SPM images
    #==========================================================================
    STAT <- xCon[[Ic[1L]]]$STAT
    #VspmSv   = cat(1,xCon(Ic).Vspm); # ONLY FOR FDR??
    VspmSv <- xCon[[Ic[1L]]]$Vspm
    stopifnot(nc == 1L)

    #-Degrees of Freedom and STAT string describing marginal distribution
    #--------------------------------------------------------------------------
    df     <- c(xCon[[Ic[1L]]]$eidf, xX$erdf)
    if(nc > 1L) {
        if(n > 1L) {
	    str <- sprintf("^{%d {Ha:k <= %d}}",nc,(nc-n)+1)
        } else {
	    str <- sprintf("^{%d {Ha:k=%d}}",nc,(nc-n)+1)
	}
    } else {
        str <- ""
    }
    
    if(STAT == "T") {
        STATstr <- sprintf("T%s_{%.0f}", str, df[2])
    } else if(STAT == "F") {
        STATstr <- sprintf("F%s_{%.0f,%.0f}", str, df[1], df[2])
    } else if(STAT == "P") {
        STATstr <- sprintf("PPM_{%0.2f}", df[1])
    }

    #-Compute conjunction as minimum of SPMs
    #--------------------------------------------------------------------------
    #for i = Ic
    #	        Z = min(Z,spm_get_data(xCon(i).Vspm,XYZ));
    #	end
    Z <- SPM$xCon[[Ic[1L]]]$Vspm[XYZ] # dimension: 1xS (vector)
    #if(length(Ic) > 1L) {
    #    for(i in 2:length(Ic)) {
    #    
    # 	}
    #}

    #-Copy of Z and XYZ before masking, for later use with FDR
    #--------------------------------------------------------------------------
    #XYZum <- SPM$XYZ
    #Zum <- Z

    #-Compute mask and eliminate masked voxels
    #--------------------------------------------------------------------------

    # TODO!!!!


    #==========================================================================
    # - H E I G H T   &   E X T E N T   T H R E S H O L D S
    #==========================================================================
    if(thresDesc == "FWE") {
        thresDesc <- paste("p<", u, " (", thresDesc, ")", sep="")
	    u <- spm_uc(u, df, STAT, R, n, S)
    } else if(thresDesc == "FDR") {
        thresDesc <- paste("p<", u, " (", thresDesc, ")", sep="")
        u <- spm_uc_FDR(u, df, STAT, n, VspmSv, 0)
    } else if(thresDesc == "none") {
        if(u < 1.0) {
            thresDesc <- paste("p<", u, " (unc.)", sep="")
	    u = spm_u(u^(1/n), df, STAT)
	} else {
	    thresDesc <- paste("STAT = ", u)
	}
    } else {
        stop("unknown correction method in getSPM: ", thresDesc)
    }

    #-Compute p-values for topological and voxel-wise FDR (all search voxels)
    #----------------------------------------------------------------------
    Ps = NA
    
    #-Peak FDR
    #----------------------------------------------------------------------
    # [up,Pp] = spm_uc_peakFDR(0.05,df,STAT,R,n,Zum,XYZum,u);
    up = NA; Pp = NA

    #-Cluster FDR
    #----------------------------------------------------------------------
    #if STAT == 'T' && n == 1
    #    V2R        = 1/prod(SPM.xVol.FWHM(SPM.xVol.DIM > 1));
    #    [uc,Pc,ue] = spm_uc_clusterFDR(0.05,df,STAT,R,n,Zum,XYZum,V2R,u);
    #else
    #    uc  = NaN;
    #    ue  = NaN;
    #    Pc  = [];
    #end
    uc = NA; ue = NA; Pc= NA

    #-Peak FWE
    #%----------------------------------------------------------------------
    uu      <- spm_uc(0.05,df,STAT,R,n,S)


    #-Calculate height threshold filtering
    #--------------------------------------------------------------------------
    Q <- which(Z > u)

    #-Apply height threshold
    #--------------------------------------------------------------------------
    Z   <- Z[Q]
    XYZ <- XYZ[Q,] 
    if(length(Q) == 0L) {
        txt <- sprintf("No voxels survive height threshold at u=%0.2g", u)
        warning(txt)
    }

    #-Extent threshold
    #--------------------------------------------------------------------------
    if(length(XYZ) > 0L && k > 0L) {
        #-Calculate extent threshold filtering
        #----------------------------------------------------------------------
        A     <- spm_clusters(XYZ)
        Q     <- integer(0)
        for(i in 1:max(A)) {
            j <- which(A == i)
            if(length(j) >= k) {
                Q <- c(Q, j)
            }
        }

        # ...eliminate voxels
        #----------------------------------------------------------------------
        Z     <- Z[Q]
        XYZ   <- XYZ[Q,]
        if(length(Q) == 0L) {
            txt <- sprintf("No voxels survive height threshold at at k=%0.2g",k)
            warning(txt)
        }
    } else {
        k <- 0
    }

    #-Assemble output structures of unfiltered data
    #=========================================================================

    xSPM <- list(
            swd       = SPM$SPMR$output.dir,
            title     = titlestr,
            Z         = Z,
            n         = n,
            STAT      = STAT,
            df        = df,
            STATstr   = STATstr,
            Ic        = Ic,
            Im        = Im,
            pm        = pm,
            Ex        = Ex,
            u         = u,
            k         = k,
            XYZ       = XYZ,
            XYZmm     = t(SPM$xVol$M[1:3,] %*% t(cbind(XYZ, 1))),
            S         = SPM$xVol$S,
            R         = SPM$xVol$R,
            FWHM      = SPM$xVol$FWHM,
            M         = SPM$xVol$M,
            iM        = SPM$xVol$iM,
            DIM       = SPM$xVol$DIM,
            VOX       = VOX,
            Vspm      = VspmSv,
            thresDesc = thresDesc)

    #-RESELS per voxel (density) if it exists
    #--------------------------------------------------------------------------
    xSPM$VRpv <- SPM$xVol$VRpv
    xSPM$units <- SPM$xVol$units

    #-p-values for topological and voxel-wise FDR
    #--------------------------------------------------------------------------
    xSPM$Ps <- Ps  # voxel   FDR
    xSPM$Pp <- Pp  # peak    FDR
    xSPM$Pc <- Pc  # cluster FDR

    #-0.05 critical thresholds for FWEp, FDRp, FWEc, FDRc
    #--------------------------------------------------------------------------
    xSPM$uc <- c(uu, up, ue, uc)
    
    xSPM
}
