#' @title spm_fMRI_design
#' @name SPM: Assembles a design for fMRI studies
#' @usage [SPM] = spm_fMRI_design(SPM)
#
# 1st level
#-------------------------------------------------------------------------
# SPM.
#
#       xY: [1x1 struct] - data structure
#    nscan: [1xs double] - nscan(s) = number of scans in session s
#      xBF: [1x1 struct] - Basis function structure
#     Sess: [1xs struct] - Session structure array
#       xX: [1x1 struct] - Design matrix structure
#
#
#    2nd level
#   ----------------------------------------------------------------------
#    SPM.xY
#           RT: - repetition time {seconds)
#
#    SPM.xBF
#            T: - number of time bins per scan
#           T0: - first time bin (see slice timing)
#        UNITS: - 'scans'|'secs' - units in which onsets are specified
#     Volterra: - 1|2 - order of [Volterra] convolution
#           dt: - length of time bin {seconds}
#         name: - name of basis set
#       length: - support of basis set {seconds}
#        order: - order of basis set
#           bf: - basis set matrix
#
#    SPM.Sess(s)
#            U: - Input structure array
#            C: - User specified covariate structure
#          row: - scan   indices for session s
#          col: - effect indices for session s
#           Fc: - F Contrast information for input-specific effects
#
#    SPM.xX
#            X: - design matrix
#           iH: - vector of H partition (indicator variables) indices
#           iC: - vector of C partition (covariates)          indices
#           iB: - vector of B partition (block effects)       indices
#           iG: - vector of G partition (nuisance variables)  indices
#         name: - cellstr of names for design matrix columns
#
#
#        3rd level
#       ------------------------------------------------------------------
#        SPM.Sess(s).U
#               dt: - time bin length {seconds}
#             name: - {1 x j} cell of names for each input or cause
#              ons: - (q x 1) onsets for q  trials {in UNITS}
#              dur: - (q x 1) durations for trials {in UNITS}
#                P: - Parameter stucture
#                u: - (t x j) inputs or stimulus function matrix
#              pst: - (1 x k) peristimulus times (seconds)
#
#
#        SPM.Sess(s).C
#
#                C: - [kx1 double] of user specified regressors
#             name: - {1xk} cellstr of regressor names
#
#
#        SPM.Sess(s).Fc
#
#                i: - F Contrast colums for input-specific effects
#             name: - F Contrast names  for input-specific effects
#
#
#            4th level
#           --------------------------------------------------------------
#            SPM.Sess(s).U(i).P(p)
#
#                 name: - parameter name
#                    P: - (q x 1) parameter matrix
#                    h: - order of polynomial expansion (0 = none)
#                    i: - sub-indices of U(i).u for plotting
#
#' @param SPM An SPM matrix.
#' @return An SPM matrix.

spm_fMRI_design <- function(SPM) {

    fMRI_T  <- SPM$xBF$T;
    fMRI_T0 <- SPM$xBF$T0;

    # separate specifications for non-replicated sessions
    # TODO: allow for rep to be TRUE?
    rep <- FALSE

    bf <- SPM$xBF$bf
    V  <- SPM$xBF$Volterra

    Xx    <- matrix(0,0,0)
    Xb    <- matrix(0,0,0)
    Xname <- ""
    Bname <- ""

    for(s in 1:length(SPM$nscan)) {
     
        # number of scans for this session
        k <- SPM$nscan[s]

        if(s == 1 || !rep) {

            # Get inputs, neuronal causes or stimulus functions U
            U <- spm_get_ons(SPM, s)

            # Convolve stimulus functions with basis functions
            X  <- spm_Volterra(U, bf, V)
            Fc <- attr(X, "Fc")
            Xn <- colnames(X)

            # Resample regressors at acquisition times (32 bin offset)
            X <- X[(0:(k - 1))*fMRI_T + fMRI_T0 + 32,,drop=FALSE]

            # and orthogonalise (within trial type)
            for(i in 1:length(Fc)) {
                X[,Fc[[i]]$i] = spm_orth( X[,Fc[[i]]$i] )
            }

            # get user specified regressors?
            C     <- SPM$Sess[[s]]$C$C;
            Cname <- SPM$Sess[[s]]$C$name;

            # append mean-corrected regressors and names
            reg_rows <- nrow(C)
            if( reg_rows > 0 ) {
                if(reg_rows != k) {
                    str1 <- "Error in spm_fMRI_design.R\n"
                    str2 <- sprintf("Session %d has %d scans but regressors have %d entries\n", s, k, reg_rows)
                    str3 <- "These numbers should match\n"
                    stop(c(str1,str2,str3))
                }
                X  <- cbind(X, scale(C, center=TRUE, scale=FALSE))
                Xn <- colnames(X)
            }

            # Confounds: Session effects
            B  <- matrix(1, nrow=k, ncol=1)
            Bn <- "constant"
        } # s == 1

        # Session structure array
        SPM$Sess[[s]]$U      = U
        SPM$Sess[[s]]$C$C    = C
        SPM$Sess[[s]]$C$name = Cname
        SPM$Sess[[s]]$row    = nrow(Xx) + (1:k)
        SPM$Sess[[s]]$col    = ncol(Xx) + (1:ncol(X))
        SPM$Sess[[s]]$Fc     = Fc

        # Append names
        for(i in 1:length(Xn)) {
            Xname <- c(Xname, sprintf("Sn(%i) ",s), Xn[i])
        }
        for(i in 1:length(Bn)) {
            Bname <- c(Bname, sprintf("Sn(%i) ",s), Bn[i])
        }

        # append into Xx and Xb
        if(s > 1) {
            ## FIXME: UNTESTED!!
            Xx <- bdiag(Xx, X)
            Xb <- bdiag(Xb, B)
        } else {
            Xx <- X
            Xb <- B
        }

    } # s

    # finished
    SPM$xX$X      = cbind(Xx, Xb)
    SPM$xX$iH     = integer(0)
    SPM$xX$iC     = ncol(Xx)
    SPM$xX$iB     = 1:ncol(Xb) + ncol(Xx)
    SPM$xX$iG     = integer(0)
    SPM$xX$name   = c(Xname, Bname)
 
    SPM
}

# R code for forming block diagonal matrix
# posted by Berten Gunter on R-help:
# Berton Gunter gunter.berton at gene.com
# Fri Sep 2 01:08:58 CEST 2005
bdiag<-function(...){
    mlist<-list(...)
    ## handle case in which list of matrices is given
    if(length(mlist)==1)mlist<-unlist(mlist,rec=FALSE)
    csdim<-rbind(c(0,0),apply(sapply(mlist,dim),1,cumsum ))
    ret<-array(0,dim=csdim[length(mlist)+1,])
    add1<-matrix(rep(1:0,2),nc=2)
    for(i in seq(along=mlist)){
        indx<-apply(csdim[i:(i+1),]+add1,2,function(x)x[1]:x[2])
          ## non-square matrix
        if(is.null(dim(indx)))ret[indx[[1]],indx[[2]]]<-mlist[[i]]
            ## square matrix
        else ret[indx[,1],indx[,2]]<-mlist[[i]]
        }
    ret
}

