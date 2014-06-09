# Fills in SPM.xCon and writes con_????.img, ess_????.img and SPM?_????.img
# FORMAT [SPM] = spm_contrasts(SPM,Ic)
#
# SPM - SPM data structure
# Ic  - indices of xCon to compute

spm_contrasts <- function(SPM, Ic=1:length(SPM$xCon)) {

    # Get contrast definitions (if available)
    xCon <- SPM$xCon

    # -Map parameter and hyperarameter files
    if(xCon[[1]]$STAT == "P") {
        # -Conditional estimators and error variance hyperparameters
        Vbeta <- SPM$VCbeta
    } else {
        #-OLS estimators and error variance estimate
        Vbeta <- SPM$Vbeta
        VHp   <- SPM$VResMS
    }

    # -Compute & store contrast parameters, contrast/ESS images, & SPM images
    for(ic in Ic) {

        # -Canonicalise contrast structure with required fields
        if( is.null(xCon[[ic]]$eidf) ) {
            X1o  <- spm_FcUtil(action="X1o", Fc=xCon[[ic]], sX=SPM$xX$xKXs)
            trMV <- spm_SpUtil("trMV", X1o, SPM$xX$V)
            trMVMV <- attr(trMV, "trMVMV"); attributes(trMV) <- NULL
            xCon[[ic]]$eidf <- trMV^2/trMVMV
        }

        # -Write contrast/ESS images?
        if( is.null(xCon[[ic]]$Vcon) ) {

            if(xCon[[ic]]$STAT == "T" || xCon[[ic]]$STAT == "P") {

                if(xCon[[ic]]$STAT == "P" &&
                   SPM$PPM$xCon[[ic]]$PSTAT == "F") {
                    # Chi^2 Bayesian inference for compound contrast
                    stop("Chi^2 Bayesian inference not implemented yet!")
                } else {
                    # -Implement contrast as sum of scaled beta images
                    Q <- which( abs(xCon[[ic]]$c) > 0 )
                    V <- Vbeta[Q]
                    # weight
                    for(j in 1:length(Q)) {
                        V[[j]] <- V[[j]] * xCon[[ic]]$c[Q[j]]
                    }
                    # sum
                    xCon[[ic]]$Vcon <- spm_add(V)
                }

            } else if(xCon[[ic]]$STAT == "F") {
                # -Implement ESS as sum of squared weighted beta images
                # -Residual (in parameter space) forming mtx
                h <- spm_FcUtil(action="Hsqr", Fc=xCon[[ic]], sX=SPM$xX$xKXs)

                # compute ESS image
                xCon[[ic]]$Vcon <- spm_resss(Vbeta, h, flags="")

            } else {
                stop("unknown STAT in xCon")
            }
        }
 
        # -Write inference SPM/PPM
        if( is.null(xCon[[ic]]$Vspm) || xCon[[ic]]$STAT == "P") {

            if(xCon[[ic]]$STAT == "T") {

                # -Compute SPM{t} image
                cB  <- xCon[[ic]]$Vcon
                l   <- VHp
                VcB <- as.numeric( t(xCon[[ic]]$c) %*% SPM$xX$Bcov %*% xCon[[ic]]$c)
                Z   <- cB/sqrt(l*VcB)

            } else if(xCon[[ic]]$STAT == "P") {
     
                # -Compute PPM{P} image
                stop("PPM not implemented yet!")

            } else if(xCon[[ic]]$STAT == "F") {

                # -Compute SPM{F} image
                MVM <- xCon[[ic]]$Vcon/trMV
                RVR <- VHp
                Z   <- MVM/RVR

            }

            # write  SPM - statistic image
            xCon[[ic]]$Vspm <- Z
        }
    } # end for i = Ic

    # place xCon back in SPM
    SPM$xCon <- xCon

    SPM
}

