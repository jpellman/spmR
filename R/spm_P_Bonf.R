#function P = spm_P_Bonf(Z,df,STAT,S,n)
# Returns the corrected P value using Bonferroni
# FORMAT P = spm_P_Bonf(Z,df,STAT,S,n)
#
# Z     - height {minium over n values}
# df    - [df{interest} df{error}]
# STAT  - Statistical field
#       'Z' - Gaussian field
#       'T' - T - field
#       'X' - Chi squared field
#       'F' - F - field
# n     - number of conjoint SPMs
# S     - Voxel count
#
# P     - corrected   P value  - P(STAT > Z)
#
#___________________________________________________________________________
#
# spm_P_Bonf returns the p-value of Z corrected by the Bonferroni
# inequality. 
#
# If n > 1 a conjunction probility over the n values of the statistic
# is returned
#
#___________________________________________________________________________
# MATLAB version: Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

spm_P_Bonf <- function(Z,df,STAT,S,n) {

   if(STAT == "Z") {
        P <- 1 - pnorm(Z)
    } else if(STAT == "T") {
        P <- 1 - pt(Z, df[2])
    } else if(STAT == "X") {
        P <- 1 - pchisq(Z, df[2])
    } else if(STAT == "F") {
        P <- 1- pf(Z, df[1], df[2])
    } else {
        stop("wrong value for STAT argument", STAT)
    }

    P <- S*P^n
    P <- min(P, 1.0)

    P
}

