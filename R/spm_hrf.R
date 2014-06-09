# returns a hemodynamic response function
# FORMAT [hrf,p] = spm_hrf(RT,[p])
# RT   - scan repeat time
# p    - parameters of the response function (two gamma functions)
#
#                                                     defaults
#                                                    (seconds)
#   p(1) - delay of response (relative to onset)         6
#   p(2) - delay of undershoot (relative to onset)      16
#   p(3) - dispersion of response                        1
#   p(4) - dispersion of undershoot                      1
#   p(5) - ratio of response to undershoot               6
#   p(6) - onset (seconds)                               0
#   p(7) - length of kernel (seconds)                   32
#
# hrf  - hemodynamic response function
# p    - parameters of the response function


spm_hrf <- function(RT, P, fMRI_T=16) {

    # default parameters
    p <- c(6, 16, 1, 1, 6, 0, 32)
    if(!missing(P)) {
        p[1:length(P)] <- P
    }

    # up-sampling
    dt    = RT/fMRI_T
    x     = 0:(p[7]/dt) - p[6]/dt

    # compute hrf
    hrf   = ( dgamma(x, shape=p[1]/p[3], rate=dt/p[3]) -
              dgamma(x, shape=p[2]/p[4], rate=dt/p[4])/p[5] )

    # down-sampling
    hrf   = hrf[0:(p[7]/RT)*fMRI_T + 1]

    # rescaling
    hrf   = as.matrix( hrf/sum(hrf) )

    # return 'p' vector as attribute 'p'
    attr(hrf, "p") <- p

    hrf
}
