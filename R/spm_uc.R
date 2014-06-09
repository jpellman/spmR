# corrected critical height threshold at a specified significance level
# FORMAT [u] = spm_uc(a,df,STAT,R,n,S)
# a     - critical probability - {alpha}
# df    - [df{interest} df{residuals}]
# STAT  - Statistical field
#       'Z' - Gaussian field
#       'T' - T - field
#       'X' - Chi squared field
#       'F' - F - field
# R     - RESEL Count {defining search volume}
# n     - number of conjoint SPMs
# S     - Voxel count
#
# u     - critical height {corrected}
#
#__________________________________________________________________________
#
# spm_uc corrected critical thresholds, using the minimum of different
# valid methods.
#
# See the individual methods for details
#
#     spm_uc_RF
#     spm_uc_Bonf
#
#__________________________________________________________________________
# MATLAB version: Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

spm_uc <- function(a,df,STAT,R,n,S=NULL) {

    u <- spm_uc_RF(a,df,STAT,R,n)

    # empty S?
    if(!is.null(S)) {
        u_Bonf <- spm_uc_Bonf(a,df,STAT,S,n)
        u <- min(c(u, u_Bonf))
    }

    u
}


# corrected critical height threshold at a specified significance level
# FORMAT [u] = spm_uc_Bonf(a,df,STAT,S,n)
# a     - critical probability - {alpha}
# df    - [df{interest} df{residuals}]
# STAT  - Statistical field
#       'Z' - Gaussian field
#       'T' - T - field
#       'X' - Chi squared field
#       'F' - F - field
#       'P' - P - value
# S     - Voxel count
# n     - number of conjoint SPMs
#
# u     - critical height {corrected}
#
#___________________________________________________________________________
#
# spm_uc returns the corrected critical threshold at a specified significance
# level (a). If n > 1 a conjunction the probability over the n values of the 
# statistic is returned.
#___________________________________________________________________________

spm_uc_Bonf <- function(a,df,STAT,S,n) {

    u <- spm_u((a/S)^(1/n),df,STAT)
    u
}

# corrected critical height threshold at a specified significance level
# FORMAT [u] = spm_uc_RF(a,df,STAT,R,n)
# a     - critical probability - {alpha}
# df    - [df{interest} df{residuals}]
# STAT  - Statistical field
#         'Z' - Gaussian field
#         'T' - T field
#         'X' - Chi-squared field
#         'F' - F field
# R     - RESEL Count {defining search volume}
# n     - number of conjoint SPMs
#
# u     - critical height {corrected}
#
#__________________________________________________________________________
#
# spm_uc returns the corrected critical threshold at a specified significance
# level (a). If n > 1 a conjunction the probability over the n values of
# the statistic is returned.
#__________________________________________________________________________

spm_uc_RF <- function(a,df,STAT,R,n) {

    # find approximate value
    #--------------------------------------------------------------------------
    u  <- spm_u((a/max(R))^(1/n),df,STAT)
    du <- 1e-6

    # approximate estimate using E{m}
    #--------------------------------------------------------------------------
    d <- 1
    while(abs(d) > 1e-6) {
        p <- spm_P_RF(1,0,u,df,STAT,R,n)[1]
        q <- spm_P_RF(1,0,u + du,df,STAT,R,n)[1] 
        d <- (a - p)/((q - p)/du)
        u <- u + d
    }

    # refined estimate using 1 - exp(-E{m})
    #--------------------------------------------------------------------------
    d  <- 1
    while(abs(d) > 1e-6) {
        p <- spm_P_RF(1,0,u,df,STAT,R,n)[1]
        q <- spm_P_RF(1,0,u + du,df,STAT,R,n)[1]
        d <- (a - p)/((q - p)/du)
        u <- u + d
    }

    u
}


#function [u] = spm_u(a,df,STAT)
# uncorrected critical height threshold at a specified significance level
# FORMAT [u] = spm_u(a,df,STAT)
# a     - critical probability - {alpha}
# df    - [df{interest} df{error}]
# STAT  - Statistical field
#               'Z' - Gaussian field
#               'T' - T - field
#               'X' - Chi squared field
#               'F' - F - field
#               'P' - P - value
#
# u     - critical height {uncorrected}
# _________________________________________________________________________
#
# spm_u returns the uncorrected critical threshold at a specified 
# significance
# _________________________________________________________________________
# Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
spm_u <- function(a,df,STAT) {

    if(STAT == "Z") {
        u <- qnorm(1-a)
    } else if(STAT == "T") {
        u <- qt(1-a, df[2])
    } else if(STAT == "X") {
        u <- qchisq(1-a, df[2])
    } else if(STAT == "F") {
        u <- qf(1-a, df[1], df[2])
    } else if(STAT == "P") {
        u <- a
    } else {
        stop("wrong value for STAT argument", STAT)
    }

    u
}

