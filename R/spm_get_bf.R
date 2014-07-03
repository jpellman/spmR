#' @name spm_get_bf
#' @title SPM: fills in basis function structure
#' @usage [xBF] = spm_get_bf(xBF);
#' @param a basis function structure
# xBF$dt      - time bin length {seconds}
# xBF$name    - description of basis functions specified
# xBF$length  - window length (secs)
# xBF$order   - order
# xBF$bf      - Matrix of basis functions
#
# xBF$name  'hrf'
#       'hrf (with time derivative)'
#       'hrf (with time and dispersion derivatives)'
#       'Fourier set'
#       'Fourier set (Hanning)'
#       'Gamma functions'
#       'Finite Impulse Response'};
#
# (any other specification will default to hrf)
#' @description spm_get_bf prompts for basis functions to model event or epoch-related
#' responses.  The basis functions returned are unitary and orthonormal
#' when defined as a function of peri-stimulus time in time-bins.
#' It is at this point that the distinction between event and epoch-related 
#' responses enters.

spm_get_bf <- function(xBF) {

    dt <- xBF$dt   
    l  <- xBF$length
    h  <- xBF$order

    # assemble basis functions
    if(xBF$name %in% c("Fourier set", 
                       "Fourier set (Hanning)",
                       "Gamma functions",
                       "Finite Impulse Response"
                       )) {
        stop(cat(xBF$name), "basis functions not supported yet!\n")
    } else if(xBF$name == "NONE") {
        bf <- 1
    } else if( grepl("hrf", xBF$name) ) {
        # canonical hemodynamic response function
        bf <- spm_hrf(dt); p <- attr(bf, "p")

        # add time derivative
        if( grepl("time", xBF$name) ) {
            dp <- 1  
            p[6] <- p[6] + dp
            D <- (bf - spm_hrf(dt,p))/dp
            bf <- cbind(bf, D)
            p[6] <- p[6] - dp

            # add dispersion derivative
            if( grepl("dispersion", xBF$name) ) {
                dp <- 0.01; p[3] <- p[3] + dp
                D <- (bf[,1] - spm_hrf(dt,p))/dp
                bf <- cbind(bf, D)
            }
        }
        xBF$length <- nrow(bf)*dt
        xBF$order  <- ncol(bf)
    }

    # Orthogonalize and fill in basis function structure
    xBF$bf <- spm_orth(bf)

    xBF
}
