# Returns the [un]corrected P value using unifed EC theory
# FORMAT [P p Ec Ek] = spm_P_RF(c,k,z,df,STAT,R,n)
#
# c     - cluster number 
# k     - extent {RESELS}
# z     - height {minimum over n values}
# df    - [df{interest} df{error}]
# STAT  - Statistical field
#       'Z' - Gaussian field
#       'T' - T - field
#       'X' - Chi squared field
#       'F' - F - field
# R     - RESEL Count {defining search volume}
# n     - number of component SPMs in conjunction
#
# P     - corrected   P value  - P(C >= c | K >= k}
# p     - uncorrected P value
# Ec    - expected number of clusters (maxima)
# Ek    - expected number of resels per cluster
#
#__________________________________________________________________________
#
# spm_P_RF returns the probability of c or more clusters with more than
# k resels in volume process of R RESELS thresholded at u.  All p values
# can be considered special cases:
#
# spm_P_RF(1,0,z,df,STAT,1,n) = uncorrected p value
# spm_P_RF(1,0,z,df,STAT,R,n) = corrected p value {based on height z)
# spm_P_RF(1,k,u,df,STAT,R,n) = corrected p value {based on extent k at u)
# spm_P_RF(c,k,u,df,STAT,R,n) = corrected p value {based on number c at k and u)
# spm_P_RF(c,0,u,df,STAT,R,n) = omnibus   p value {based on number c at u)
#
# If n > 1 a conjunction probility over the n values of the statistic
# is returned
#__________________________________________________________________________
#
# References:
# 
# [1] Hasofer AM (1978) Upcrossings of random fields
# Suppl Adv Appl Prob 10:14-21
# [2] Friston KJ et al (1994) Assessing the Significance of Focal Activations
# Using Their Spatial Extent
# Human Brain Mapping 1:210-220
# [3] Worsley KJ et al (1996) A Unified Statistical Approach for Determining
# Significant Signals in Images of Cerebral Activation
# Human Brain Mapping 4:58-73
#__________________________________________________________________________
# MATLAB version: Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

# utility function taken from Biodem package
# take the 'n-th' power of matrix X
mx.exp <- function (X, n) {
    if (n != round(n)) {
        n <- round(n)
        warning("rounding exponent `n' to", n)
    }
    phi <- diag(nrow = nrow(X))
    pot <- X
    while (n > 0) {
        if (n%%2) 
            phi <- phi %*% pot
        n <- n%/%2
        pot <- pot %*% pot
    }
    return(phi)
}


spm_P_RF <- function(c,k,Z,df,STAT,R,n) {

    eps <- .Machine$double.eps

    # get expectations
    #==========================================================================

    # get EC densities
    #--------------------------------------------------------------------------
    D <- which(R > 0)[length(R)]
    R <- R[1:D]
    G <- sqrt(pi)/gamma((1:D)/2)
    EC <- spm_ECdensity(STAT,Z,df)
    EC <- pmax(EC[1:D], eps)

    # corrected p value
    #--------------------------------------------------------------------------
    P   <- toeplitz( as.numeric(t(EC) * G) ); P[lower.tri(P)] <- 0.0
    P   <- mx.exp(P, n)
    P   <- P[1,]
    EM  <- (R/G)*P         # <maxima> over D dimensions
    Ec  <- sum(EM)         # <maxima>
    EN  <- P[1]*R[D]       # <resels>
    Ek  <- EN/EM[D]        # Ek = EN/EM(D);

    # get P{n > k}
    #==========================================================================

    # assume a Gaussian form for P{n > k} ~ exp(-beta*k^(2/D))
    # Appropriate for SPM{Z} and high d.f. SPM{T}
    #--------------------------------------------------------------------------

    D <- D - 1
    if(k == 0 || D == 0) {
        p <- 1
    } else if(STAT == "Z") {
        beta <- (gamma(D/2 + 1)/Ek)^(2/D)
        p <- exp(-beta*(k^(2/D)))
    } else if (STAT=='T'){
        beta <- (gamma(D/2 + 1)/Ek)^(2/D)
        p <- exp(-beta*(k^(2/D)))
    } else if (STAT=='X'){
        beta <- (gamma(D/2 + 1)/Ek)^(2/D)
        p <- exp(-beta*(k^(2/D)))
    } else if (STAT=='F'){
        beta <- (gamma(D/2 + 1)/Ek)^(2/D)
        p <- exp(-beta*(k^(2/D)))
    } else {
        cat("STAT incorrectly specified")
    }

    # Poisson clumping heuristic {for multiple clusters}
    #==========================================================================
    P <- 1 - ppois(c - 1,(Ec + eps)*p);

    # set P and p = [] for non-implemented cases
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(k > 0 && (STAT == "X" || STAT == "F")) {
        P <- a.numeric(0); p <- as.numeric(0)
    }
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    out <- c(P, p, Ec, Ek)
    out
} 

#==========================================================================
# spm_ECdensity
#==========================================================================
#unction [EC] = spm_ECdensity(STAT,t,df)
# Returns the EC density
#__________________________________________________________________________
#
# Reference : Worsley KJ et al 1996, Hum Brain Mapp. 4:58-73
#
#--------------------------------------------------------------------------
spm_ECdensity <- function(STAT,t,df) {

    # EC densities (EC)
    #--------------------------------------------------------------------------
    #t = t(:)';

    EC <- matrix(rep(NA,4),nrow=4)

    if(STAT=="Z") {
        # Gaussian Field
        #----------------------------------------------------------------------
        a <- 4*log(2)
        b <- exp(-t^2/2)

        EC[1,] <- 1 - pnorm(t)
        EC[2,] <- a^(1/2)/(2*pi)*b
        EC[3,] <- a/((2*pi)^(3/2))*b*t
        EC[4,] <- a^(3/2)/((2*pi)^2)*b*(t^2-1)

    } else if (STAT=="T") {
        # T - Field
        #----------------------------------------------------------------------
        v <- df[2]
        a <- 4*log(2)
        b <- exp(lgamma((v+1)/2) - lgamma(v/2))
        c <- (1+t^2/v)^((1-v)/2)
  
        EC[1,] <- 1 - pt(t,v)
        EC[2,] <- a^(1/2)/(2*pi)*c
        EC[3,] <- a/((2*pi)^(3/2))*c*t/((v/2)^(1/2))*b
        EC[4,] <- a^(3/2)/((2*pi)^2)*c*((v-1)*(t^2)/v-1)

     } else if (STAT=="X") {
        # X - Field
        #----------------------------------------------------------------------
        v <- df[2]
        a <- (4*log(2))/(2*pi)
        b <- t^(1/2*(v-1))*exp(-t/2-lgamma(v/2))/2^((v-2)/2)

        EC[1,] <- 1 - pchisq(t,v)
        EC[2,] <- a^(1/2)*b
        EC[3,] <- a*b*(t-(v-1))
        EC[4,] <- a^(3/2)*b*(t^2-(2*v-1)*t+(v-1)*(v-2))

    } else if (STAT=="F") {
        # F Field
        #----------------------------------------------------------------------
        k <- df[1] ### CHECKME!
        v <- df[2] ### CHECKME!
        a <- (f*log(2))/(2*pi)
        b <- lgamma(v/2) + lgamma(k/2)

        EC[1,] <- 1 - pf(t, df[1], df[2])
        EC[2,] <- a^(1/2)*exp(lgamma((v+k-1)/2)-b)*2^(1/2)*(k*t/v)^(1/2*(k-1))*(1+k*t/v)^(-1/2*(v+k-2))
        EC[3,] <- a*exp(lgamma((v+k-2)/2)-b)*(k*t/v)^(1/2*(k-2))*(1+k*t/v)^(-1/2*(v+k-2))*((v-1)*k*t/v-(k-1))
        EC[4,] <- a^(3/2)*exp(lgamma(((v+k-3)/2))-b)*2^(-1/2)*(k*t/v)^(1/2*(k-3))*(1+k*t/v)^(-1/2*(v+k-2))*((v-1)*(v-2)*(k*t/v)^2-2(2*v*k-v-k-1)*(k*t/v)+(k-1)*(k-2))
    } else {
        stop("wrong value for argument STAT:", STAT)
    }

    EC
}
