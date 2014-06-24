#' @name spm_P
#function [P,p,Ec,Ek] = spm_P(c,k,Z,df,STAT,R,n,S)
#' @title SPM: Returns the [un]corrected P value using unifed EC theory
#' @usage FORMAT [P p Ec Ek] = spm_P(c,k,Z,df,STAT,R,n,S)
#
#' @param c - cluster number 
#' @param k - extent {RESELS}
#' @param Z - height {minimum over n values}
#' @param df - [df{interest} df{error}]
#' @param STAT  - Statistical field: 'Z' - Gaussian field, 'T' - T - field, 'X' - Chi squared field, 'F' - F - field,  'P' - Posterior probability
#' @param R - RESEL Count {defining search volume}
#' @param n - number of component SPMs in conjunction
#' @param S - Voxel count
#
#' @return A list containing: P     - corrected   P value - P(C >= c | K >= k}, p     - uncorrected P value, 
#' Ec    - expected total number of clusters and Ek    - expected total number of resels per cluster
#
#__________________________________________________________________________
#
#' @description spm_P determines corrected and uncorrected p values, using the minimum
#' of different valid methods. 
#'
#' See the individual methods for details
#
#' @seealso spm_P_RF
#' @seealso spm_P_Bonf
#
#__________________________________________________________________________
# MATLAB version: Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

spm_P <- function(c,k,Z,df,STAT,R,n,S=NULL) {

    out <- spm_P_RF(c,k,Z,df,STAT,R,n)
    P <- out[1]
    p <- out[2]
    Ec <- out[3]
    Ek <- out[4]

    # Use lower Bonferroni P value (if possible)
    #==========================================================================
    if(!is.null(S) && (c == 1 && k == 0) && !(length(R) == 1L && R == 1)) {
        P = min(P, spm_P_Bonf(Z,df,STAT,S,n))
    }

    c(P, p, Ec, Ek)
}

