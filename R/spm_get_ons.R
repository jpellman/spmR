#' @name spm_get_ons
#' @title SPM: returns input [designed effects] structures
#' @usage [U] = spm_get_ons(SPM,s)
#
#' @param SPM The SPM object.
#' @param s session number (used by batch system)
#
#' @return U (1 x n)   struct array of (n) trial-specific structures
#
#' @description   U(i).name   - cell of names for each input or cause
#'   U(i).u      - inputs or stimulus function matrix
#'   U(i).dt     - time bin (seconds)
#'   U(i).ons    - onsets    (in SPM.xBF.UNITS)
#'   U(i).dur    - durations (in SPM.xBF.UNITS)
#'   U(i).P      - parameter struct.
#'
#'       U(i).P(p).name - parameter name
#'       U(i).P(p).P    - parameter vector
#'       U(i).P(p).h    - order of polynomial expansion
#'       U(i).P(p).i    - sub-indices of u pertaining to P

spm_get_ons <- function(SPM, s=1) {

    k     <- SPM$nscan[s]
    T     <- SPM$xBF$T
    dt    <- SPM$xBF$dt
    UNITS <- SPM$xBF$UNITS

    if(UNITS == "scans") {
        TR <- T*dt
    } else if(UNITS == "secs") {
        TR <- 1
    } else {
        stop("Unkown value in UNITS")
    }

    # get inputs and names
    U <- SPM$Sess[[s]]$U
    v <- length(U)

    # get trials
    for(i in 1:v) {
        # get names
        Uname <- U[[i]]$name[1]

        # 1. get main [trial] effect

        # 1.a onsets
        ons <- U[[i]]$ons

        # 1.b durations
        dur <- U[[i]]$dur

        # 1.c peri-stimulus times {seconds}
        pst <- 0:(k-1)*T*dt - min(ons)*TR
        for(j in 1:length(ons)) {
            w <- 0:(k-1)*T*dt - ons[j]*TR
            v.idx <- which(w >= 0)
            pst[v.idx] <- w[v.idx]
        }

        # 2. add parameters x trial interactions

        # 2.a get parameter stucture xP
        xP <- U[[i]]$P

        # 2.b interaction with causes (u) - 1st = main effects
        u <- as.matrix(ons^0)
        if(length(xP) > 0) {
            for(q in 1:length(xP)) {
                xP[[q]]$i = c(1, 1:xP[[q]]$h + ncol(u))
                for(j in 1:xP[[q]]$h) {
                    u   = c(u, xP[[q]]$P.^j)
                    str = sprintf("%sx%s^%d", Uname, xP[[q]]$name, j)
                    Uname <- c(Uname, str)
                }
            }
        }
 
        # orthogonalize inputs
        u <- spm_orth(u)

        # and scale so sum(u*dt) = number of events, if event-related
        if(sum(dur) == 0) {
            u <- u/dt
        }

      
        # 3. create stimulus functions (32 bin offset)
        ton  <- floor(ons*TR/dt + 0.5) + 33         # onsets
        tof  <- floor(dur*TR/dt + 0.5) + ton + 1    # offset
        sf   <- matrix(0, nrow=(k*T + 128), ncol=ncol(u))
        ton  <- pmax(ton,1)
        tof  <- pmax(tof,1)
        for(j in 1:length(ton)) {
            if(length(sf) > ton[j]) {
                sf[ton[j],] = sf[ton[j],] + u[j,]
            }
            if(length(sf)>tof[j]) {
                sf[tof[j],] = sf[tof[j],] - u[j,]
            }
        }
        sf <- apply(sf, 2, cumsum)   # integrate
        sf <- sf[1:(k*T + 32),] # stimulus


        # fill in structure
        U[[i]]$name <- Uname         # - input names
        U[[i]]$dt   <- dt            # - time bin {seconds}
        U[[i]]$u    <- sf            # - stimulus function matrix
        U[[i]]$pst  <- pst           # - pst (seconds)
        U[[i]]$P    <- xP            # - parameter struct

    } # v

    U
}

