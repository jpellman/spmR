# This file is based on the SPM software (http://www.fil.ion.ucl.ac.uk/spm/)
# Original Matlab code: Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
# R port by Yves Rosseel (http://www.da.ugent.be)
# Please note: the R port is only an approximation; therefore, the results may
# not be identical to the SPM results; moreover, if you find any bugs in the
# R code, please blame us (http://www.da.ugent.be), not the FIL group

# Removes low frequency confounds X0
# FORMAT [Y] = spm_filter(K,Y)
# FORMAT [K] = spm_filter(K)
#
# K           - filter matrix or:
# K(s)        - struct array containing partition-specific specifications
#
# K(s).RT     - observation interval in seconds
# K(s).row    - row of Y constituting block/partition s
# K(s).HParam - cut-off period in seconds
#
# K(s).X0     - low frequencies to be removed (DCT)
# 
# Y           - data matrix
#
# K           - filter structure
# Y           - filtered data
#____________________________________________________________________
#
# spm_filter implements high-pass filtering in an efficient way by
# using the residual forming matrix of X0 - low frequency confounds
#.spm_filter also configures the filter structure in accord with the 
# specification fields if called with one argument
#____________________________________________________________


spm_filter <- function(K, Y=NULL) {

    if(is.null(Y)) {
        # set K$X0
        for(s in 1:length(K)) {
            # make high pass filter
            k  <- length(K[[s]]$row)
            n  <- floor( 2*(k*K[[s]]$RT)/K[[s]]$HParam + 1)
            X0 <- spm_dctmtx(k, n)
            K[[s]]$X0 <- X0[,-1]
        } 
        return(K)

    } else {
        if(is.list(K)) {
            for(s in 1:length(K)) {
                # select data
                y <- Y[ K[[s]]$row, ]
                
                # apply high pass filter
                y <- y - K[[s]]$X0 %*% ( t(K[[s]]$X0) %*% y )

                # reset filtered data in Y
                Y[ K[[s]]$row, ] <- y
            }
        } else { # K is simply a filter matrix
            Y <- K %*% Y
        }
        
        return(Y)
    }
}
