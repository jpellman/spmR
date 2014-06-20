#' @name spm_list
#' @title spm_list
#
#' @param action
#' @param xSPM
#' @param Num
#' @param Dis
#
#
# TabDat - Structure containing table data
#        - fields are
# .tit   - table Title (string)
# .hdr   - table header (2x12 cell array)
# .fmt   - fprintf format strings for table data (1x12 cell array)
# .str   - table filtering note (string)
# .ftr   - table footnote information (5x2 cell array)
# .dat   - table data (Nx12 cell array)

#' @description
#' spm_list characterizes SPMs (thresholded at u and k) in terms of
#' excursion sets (a collection of face, edge and vertex connected
#' subsets or clusters).  The corrected significance of the results are
#' based on set, cluster and voxel-level inferences using distributional
#' approximations from the Theory of Gaussian Fields.  These
#' distributions assume that the SPM is a reasonable lattice
#' approximation of a continuous random field with known component field
#' smoothness.
#'
#' The p values are based on the probability of obtaining c, or more,
#' clusters of k, or more, resels above u, in the volume S analysed =
#' P(u,k,c).  For specified thresholds u, k, the set-level inference is
#' based on the observed number of clusters C, = P(u,k,C).  For each
#' cluster of size K the cluster-level inference is based on P(u,K,1)
#' and for each voxel (or selected maxima) of height U, in that cluster,
#' the voxel-level inference is based on P(U,0,1).  All three levels of
#' inference are supported with a tabular presentation of the p values
#' and the underlying statistic:
#'
#' Set-level     - c    = number of suprathreshold clusters
#'               - P    = prob(c or more clusters in the search volume)
#'
#' Cluster-level - k    = number of voxels in this cluster
#'               - Pc   = prob(k or more voxels in the search volume)
#'               - Pu   = prob(k or more voxels in a cluster)
#'               - Qc   = lowest FDR bound for which this cluster would be
#'                        declared positive
#'
#' Peak-level    - T/F  = Statistic upon which the SPM is based
#'               - Ze   = The equivalent Z score - prob(Z > Ze) = prob(t > T)
#'               - Pc   = prob(Ze or higher in the search volume)
#'               - Qp   = lowest FDR bound for which this peak would be
#'                        declared positive
#'               - Pu   = prob(Ze or higher at that voxel)
#'
#' Voxel-level   - Qu   = Expd(Prop of false positives among voxels >= Ze)
#'
#' x,y,z (mm)    - Coordinates of the voxel
#'
#' The table is grouped by regions and sorted on the Ze-variate of the
#' primary maxima.  Ze-variates (based on the uncorrected p value) are the
#' Z score equivalent of the statistic. Volumes are expressed in voxels.
#'
#' Clicking on values in the table returns the value to the Matlab
#' workspace. In addition, clicking on the co-ordinates jumps the
#' results section cursor to that location. The table has a context menu
#' (obtained by right-clicking in the background of the table),
#' providing options to print the current table as a text table, or to
#' extract the table data to the Matlab workspace.
#'
#' @seealso spm_getSPM
# (see spm_getSPM for further details of xSPM structures)
#' @usage
#  Display and analysis of SPM{.}: 
#  function varargout = spm_list(varargin)
#  
#' Summary list of local maxima for entire volume of interest:
#' TabDat = spm_list('List',SPM,hReg,[Num,Dis,Str])
#'
#' List of local maxima for a single suprathreshold cluster:
#' TabDat = spm_list('ListCluster',SPM,hReg,[Num,Dis,Str])

spm_list <- function(action="List", xSPM=NULL, Num=3L, Dis=8) {

    stopifnot(action == "List")

    #-Tolerance for p-value underflow, when computing equivalent Z's
    #----------------------------------------------------------------------
    tol = .Machine$double.eps*10

    #-Extract data from xSPM
    #----------------------------------------------------------------------
    S         <- xSPM$S
    VOX       <- xSPM$VOX
    DIM       <- xSPM$DIM
    n         <- xSPM$n
    STAT      <- xSPM$STAT
    df        <- xSPM$df
    u         <- xSPM$u
    M         <- xSPM$M
    k         <- xSPM$k
    QPs       <- xSPM$Ps
    QPp       <- xSPM$Pp
    QPc       <- xSPM$Pc
    thresDesc <- xSPM$thresDesc
    
    if(STAT != "P") {
        R         <- xSPM$R
        FWHM      <- xSPM$FWHM
    }

    if(!is.null(xSPM$units)) {
        units <- xSPM$units
    } else {
        units <- c("mm", "mm", "mm")
    }

    DIM       <- DIM > 1              # non-empty dimensions
    D         <- sum(DIM)             # highest dimension
    VOX       <- VOX[DIM]             # scaling

    if(STAT != "P") {
        FWHM  <- FWHM[DIM]            # Full width at max/2
        FWmm  <- FWHM*VOX             # FWHM {units}
        V2R   <- 1/prod(FWHM)         # voxels to resels
        k     <- k*V2R                # extent threshold in resels
        R     <- R[1:(D + 1)]         # eliminate null resel counts
        if(!is.null(QPs)) QPs <- sort(QPs)  # Needed for voxel   FDR
        if(!is.null(QPp)) QPp <- sort(QPp)  # Needed for peak    FDR
        if(!is.null(QPc)) QPc <- sort(QPc)  # Needed for cluster FDR
    }

    #-Get number and separation for maxima to be reported
    #----------------------------------------------------------------------
    if(STAT == "P") {
        Title <- "Posterior Probabilities"
    } else {
        Title  <- "p-values adjusted for search volume"
    }

    #-Table header & footer
    #======================================================================
    Title <- paste("Statistics: ", Title, sep="")

    
    #-Volume, resels and smoothness (if classical inference)
    #----------------------------------------------------------------------
    TabDat.ftr <- character(9)
    if(STAT != "P") {
        #------------------------------------------------------------------
        Pz  <- spm_P(1,0,u,df,STAT,1,n,S)[1]
        Pu  <- spm_P(1,0,u,df,STAT,R,n,S)[1]
        out <- spm_P(1,k,u,df,STAT,R,n,S)
        P <- out[1]; Pn <- out[2]; Ec <- out[3]; Ek <- out[4]

        TabDat.ftr[1] <-
            sprintf("Height threshold: %s = %0.2f, p = %0.3f (%0.3f)",
            STAT,u,Pz,Pu)
        TabDat.ftr[2] <-
            sprintf("Extent threshold: k = %0.0f voxels, p = %0.3f (%0.3f)",
            k/V2R,Pn,P)
        TabDat.ftr[3] <-
            sprintf("Expected voxels per cluster, <k> = %0.3f",Ek/V2R)
        TabDat.ftr[4] <-
            sprintf("Expected number of clusters, <c> = %0.2f",Ec*Pn)
        TabDat.ftr[5] <-
            sprintf("FWEp: %0.3f, FDRp: %0.3f, FWEc: %0.0f, FDRc: %0.0f",
                    xSPM$uc[1], xSPM$uc[2], xSPM$uc[3], xSPM$uc[4])
        TabDat.ftr[6] <-
            sprintf("Degrees of freedom = [%0.1f, %0.1f]",df[1], df[2])
        TabDat.ftr[7] <- paste("FWHM = ", 
            sprintf("%0.1f %0.1f %0.1f %s %s %s;", FWmm[1], FWmm[2], FWmm[3],
                    units[1], units[2], units[3]),
            " voxels = ", 
            sprintf("%0.1f %0.1f %0.1f {voxels}", FWHM[1], FWHM[2], FWHM[3]),
            sep="")
        TabDat.ftr[8] <-
            sprintf("Volume: %0.0f = %0.0f voxels = %0.1f resels", 
            S*prod(VOX),S,R[length(R)])
        TabDat.ftr[9] <- paste("Voxel size: ", 
                             sprintf("%0.1f %0.1f %0.1f %s %s %s;", 
                                     VOX[1], VOX[2], VOX[3],
                                     units[1], units[2], units[3]),
                             sprintf("(resel = %0.2f voxels)",prod(FWHM)),
                             sep="")
    } 


    #-Characterize excursion set in terms of maxima
    # (sorted on Z values and grouped by regions)
    #======================================================================
    
    #-Workaround in spm_max for conjunctions with negative thresholds
    #----------------------------------------------------------------------
    minz      <- abs(min(xSPM$Z))
    zscores   <- 1.0 + minz + xSPM$Z
    out       <- spm_max(zscores, xSPM$XYZ)
    N <- out$N; Z <- out$Z; XYZ <- out$M; A <- out$A; L <- out$XYZ
    Z         <- Z - minz - 1.0

    #-Convert cluster sizes from voxels (N) to resels (K)
    #----------------------------------------------------------------------
    c      <- max(A)                                   #-Number of clusters
     
    # TODO
    #
    K <- N*V2R

    #-Convert maxima locations from voxels to mm
    #----------------------------------------------------------------------
    XYZmm <- t(M[1:3,] %*% t(cbind(XYZ, 1)))


    #-Table proper (& note all data in cell array)
    #======================================================================

    #-Set-level p values {c} - do not display if reporting a single cluster
    #----------------------------------------------------------------------
    if(STAT != "P") {
        Pc <- spm_P(c,k,u,df,STAT,R,n,S)[1]            # -Set-level p-value
    } else {
        Pc <- as.numeric(NA)
    }

    # we will store everything in a matrix, converting later to dataframe
    # we fill the matrix row by row
    TABLE <- matrix(as.numeric(NA), nrow=0, ncol=14L)

    #-Local maxima p-values & statistics
    #----------------------------------------------------------------------
    while(sum(is.finite(Z)) > 0L) {
        #-Find largest remaining local maximum
        #------------------------------------------------------------------
        U <- max(Z, na.rm=TRUE); i <- which(Z == U)   # largest maxima
        j <- which(A == A[i])                         # maxima in cluster

        #-Compute cluster {k} and peak-level {u} p values for this cluster
        #------------------------------------------------------------------
        if(STAT != "P") {
           
            # p-values (FWE)
            #--------------------------------------------------------------
            Pz      <- spm_P(1,0,   U,df,STAT,1,n,S)[1]  # uncorrected p value
            Pu      <- spm_P(1,0,   U,df,STAT,R,n,S)[1]  # FWE-corrected {based on Z}
            out <- spm_P(1,K[i],u,df,STAT,R,n,S)   # [un]corrected {based on K}
            Pk <- out[1]; Pn <- out[2]
           

            # q-values (FDR)
            #--------------------------------------------------------------
            #if topoFDR
            #    Qc  = spm_P_clusterFDR(K(i),df,STAT,R,n,u,QPc); % based on K
            #    Qp  = spm_P_peakFDR(U,df,STAT,R,n,u,QPp);       % based on Z
            #    Qu  = [];
            #else
            #    Qu  = spm_P_FDR(U,df,STAT,n,QPs);     % voxel FDR-corrected
            #    Qc  = [];
            #    Qp  = [];
            #end
            Qu <- Qc <- Qp <- as.numeric(NA)

            # Equivalent Z-variate
            #--------------------------------------------------------------
            if( Pz < tol) {
                Ze  <- Inf
            } else {
                Ze  <- qnorm(1 - Pz)
            }
        } else {
            Pz <- Pu <- Qu <- Pk <- Pn <- Qc <- Qp <- as.numeric(NA)
            Ze <- qnorm(U)
        }
  
        # fill in row
        TABLE <- rbind(TABLE,
        c(NA, NA, Pk, Qc, N[i], Pn, Pu, Qp, U, Ze, Pz, XYZmm[i,]))

        #-Print Num secondary maxima (> Dis mm apart)
        #------------------------------------------------------------------
        out <-  sort.int(-Z[j], index.return=TRUE)
        l <- out$x; q <- out$ix
        D <- i
        for(i in 1:length(q)) {
            d <- j[q[i]]
            DIST <- sqrt(apply(scale(XYZmm[D,,drop=FALSE], center=XYZmm[d,], scale=FALSE)[,,drop=FALSE]^2, 1, sum))
            if( min(DIST) > Dis ) {
                if(length(D) < Num) {
                    # voxel-level p values {Z}
                    #------------------------------------------------------
                    if(STAT != "P") {
                        Pz    = spm_P(1,0,Z[d],df,STAT,1,n,S)[1]
                        Pu    = spm_P(1,0,Z[d],df,STAT,R,n,S)[1]
                        #if topoFDR
                        #    Qp = spm_P_peakFDR(Z(d),df,STAT,R,n,u,QPp);
                        #    Qu = [];
                        #else
                        #    Qu = spm_P_FDR(Z(d),df,STAT,n,QPs);
                        #    Qp = [];
                        #end
                        Qu <- Qp <- as.numeric(NA)
                        if(Pz < tol) {
                            Ze <- Inf
                        } else {
                            Ze <- qnorm(1 - Pz)
                        }
                    } else {
                        Pz <- Pu <- Qu <- Qp <- as.numeric(NA)
                        Ze <- qnorm(Z[d])
                    }
                    D <- c(D, d)

                    # fill in row
                    TABLE <- rbind(TABLE,
                    c(NA, NA, NA, NA, NA, NA, Pu, Qp, Z[d], Ze, Pz, XYZmm[i,]))
                }
            }
        }
        Z[j] <- as.numeric(NA)
    }

    TABLE[1,c(1,2)] <- c(Pc, c)
    TABLE <- as.data.frame(TABLE)
    names(TABLE) <- c("p","c",
                      "cl.p.FWE.corr",
                      "cl.q.FDR.corr",
                      "cl.kE",
                      "cl.p.unc",
                      "p.FWE.corr",
                      "q.FDR.corr",
                      "T",
                      "(Z)",
                      "p.unc",
                      "X", "Y", "Z")

               
    attr(TABLE, "footer") <- TabDat.ftr                       

    TABLE
}
