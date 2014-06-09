spmr_read_nifti <- function(filename="", verbose=TRUE) {

    # read in images
    require(Rniftilib)
    if(verbose) cat("Reading in 4D nifti file... ")
    Y <- nifti.image.read(filename)
    if(verbose) {
        cat("done.\n")
        cat("Dimensions: [", dim(Y), "]\n", sep=" ")
    }

    # repair slope
    if(Y$scl.slope == 0) {
        Y$scl.slope <- 1
    }

    #if(verbose) cat("Reading in 4D nifti file... ")
    #require(oro.nifti, quietly=TRUE)
    #VY <- readNIfTI(filename)
    #Y <- VY@.Data * VY@scl_slope + VY@scl_inter
    #if(verbose) {
    #    cat("done.\n")
    #    cat("Dimensions: [", dim(Y), "]\n", sep=" ")
    #}

    Y
}


spmr_write_nifti <- function(V, output.file, check=1, datatype=16) {

    require(Rniftilib)

    # check input
    if(!is.array(V)) {
        stop("V is not an array")
    }
    if(missing(output.file)) {
        stop("please specify filename for writing nifti image")
    } else if(!is.character(output.file)) {
        stop("output.file is not a character")
    }

    # create new nifti object
    IMG <- nifti.image.new()

    # change dimensions
    IMG$dim <- dim(V)

    # change filename
    check <- nifti.set.filenames(nim=IMG, prefix=output.file,
                                 check=check,
                                 set_byte_order=1) ## FIXME

    # change datatype  - default = 16 float (32 bits/voxel)
    nifti.image.setdatatype(nim=IMG, value=datatype)

    # assign data
    if(length(dim(V) == 1)) {
        IMG[] <- V
    } else if(length(dim(V) == 2)) {
        IMG[,] <- V
    } else if(length(dim(V) == 3)) {
        IMG[,,] <- V
    } else if(length(dim(V) == 4)) {
        IMG[,,,] <- V
    } else if(length(dim(V) == 5)) {
        IMG[,,,,] <- V
    } else if(length(dim(V) == 6)) {
        IMG[,,,,,] <- V
    }

    out <- nifti.image.write(IMG)

    invisible(out)
}


spmr_fmri_specify_1stlevel <- function(RT,
                                       units,
                                       fmri_t, fmri_t0,
                                       bf.name,
                                       bf.length=NA, bf.order=NA,
                                       bf.volt,
                                       sess=list(),
                                       fact=NULL,
                                       global,
                                       cvi,
                                       mask=NULL,
                                       output.dir="/tmp",
                                       verbose=TRUE
                                      ) {


    # create empty SPM object
    SPM <- list(xY=list(RT=0,P=character(),VY=NULL), 
                nscan=numeric( length(sess) ), 
                xBF=list(), 
                Sess=list(),
                xX=list(K=list()), 
                xGX=list(iGXcalc="",sGXcalc="",
                         rg=numeric(),GM=NA, gSF<-numeric()), 
                xVi=list(form="", V=NULL, Vi=NULL), 
                xM=list(T=numeric(), TH=numeric(), I=0, VM <- list(),
                        xs <- list()), 
                xsDes=list(),
                SPMR=list(output.dir=output.dir, verbose=verbose)
               )

    # create output dir -- we will override everything!!
    dir.create(SPM$SPMR$output.dir, showWarnings=FALSE)

    # check input + insert default values
    if(missing(RT)) stop("RT value is missing")
    if(missing(units)) stop("units value is missing")
    if(missing(fmri_t))  fmri_t  <- defaults.stats.fmri.fmri_t
    if(missing(fmri_t0)) fmri_t0 <- defaults.stats.fmri.fmri_t0
    if(missing(bf.name)) {
        if(defaults.stats.fmri.hrf.derivs[1] == 0 &&
           defaults.stats.fmri.hrf.derivs[2] == 0) { 
            bf.name <- "hrf"
        } else if(defaults.stats.fmri.hrf.derivs[1] == 1 &&
                  defaults.stats.fmri.hrf.derivs[2] == 0) { 
            bf.name <- "hrf (with time derivative)"
        } else if(defaults.stats.fmri.hrf.derivs[1] == 1 &&
                  defaults.stats.fmri.hrf.derivs[2] == 1) {
            bf.name <- "hrf (with time and dispersion derivatives)"
        } else {
            stop("inconsistent values for defaults.stats.fmri.hrf.derivs")
        } 
    }
    if(missing(bf.volt)) bf.volt <- defaults.stats.fmri.volt
    if(missing(global)) global <- defaults.stats.fmri.global
    if(missing(cvi)) cvi <- defaults.stats.fmri.cvi
    if(!is.list(sess) || length(sess) == 0) stop("wrong sess argument")

    # insert default values in session/cond (eg. hpf tmod pmod)
    for(i in 1:length(sess)) {
        if( is.null(sess[[i]]$hpf) ) {
            sess[[i]]$hpf <- defaults.stats.fmri.hpf
        }
        if(!is.list(sess[[i]]$cond) || 
            length(sess[[i]]$cond) == 0) stop("wrong cond argument")
        for(j in 1:length(sess[[i]]$cond)) {
            if( is.null(sess[[i]]$cond[[j]]$tmod) ) {
                sess[[i]]$cond[[j]]$tmod <- defaults.stats.fmri.cond.tmod
            }
            if( is.null(sess[[i]]$cond[[j]]$pmod) ) {
                sess[[i]]$cond[[j]]$pmod <- list()
            }
        }
        if( is.null(sess[[i]]$regress) ) {
            sess[[i]]$regress <- list()
        }
    }


    ## we first read in the 4D nifti files: one per session/subject
    ###########
    ## FIXME!!!
    ## how to 'glue' 4D images into one, along time dimension?
    #SPM$SPMR$VY <- VY[[1]] 
    ###########

    # this is just a pointer to the data, per session
    SPM$SPMR$VY <- vector("list", length=length(sess))
    for(i in 1:length(sess)) {
        SPM$SPMR$VY[[i]] <- spmr_read_nifti(sess[[i]]$scans)
    }

    ##
    ## FILE: SPMROOT/config/spm_run_fmri_spec.m 
    ##
    SPM$xY$RT <- RT

    # basis function variables
    xBF <- list()
    xBF$UNITS    <- units
    xBF$T        <- fmri_t
    xBF$T0       <- fmri_t0
    xBF$dt       <- RT/fmri_t
    xBF$name     <- bf.name
    xBF$length   <- bf.length
    xBF$order    <- bf.order
    xBF$Volterra <- bf.volt
    xBF          <- spm_get_bf(xBF)

    SPM$xBF <- xBF

    # session info (TODO: multicondition/multiregression from file!)
    for(i in 1:length(sess)) {
        ## FIXME: nscans is function of number of files!
        ## SPM$nscan[i]  <- length( sess[[i]]$scans )
        SPM$nscan[i]  <- sess[[i]]$nscans

        SPM$xY$P      <- c(SPM$xY$P, sess[[i]]$scans)
        SPM$Sess[[i]] <- list()
 
        # condition info
        U <- vector("list", length=length(sess[[i]]$cond))
        for(j in 1:length(sess[[i]]$cond)) {
            U[[j]] <- list()
            U[[j]]$name <- sess[[i]]$cond[[j]]$name
            U[[j]]$ons  <- sess[[i]]$cond[[j]]$onset
            U[[j]]$dur  <- sess[[i]]$cond[[j]]$duration
            if(length(U[[j]]$dur) == 1) {
                U[[j]]$dur <- U[[j]]$dur * rep(1, length(U[[j]]$ons))
            } else if(length(U[[j]]$dur) != length(U[[j]]$ons)) {
                stop("Mismatch between number of onset and number of durations")
            }

            P <- list()
            q1 <- 0
            # time effects? 
            if(sess[[i]]$cond[[j]]$tmod > 0) {
                P[[1]] <- list()
                P[[1]]$name <- "time"
                P[[1]]$P    <- U[[j]]$ons * RT/60
                P[[1]]$h    <- sess[[i]]$cond[[j]]$tmod
                q1          <- 1
            }
  
            # parametric modulation?
            if( length(sess[[i]]$cond[[j]]$pmod) > 0 ) {
                for(q in 1:length(sess[[i]]$cond[[j]]$pmod)) {
                    q1 <- q1 + 1
                    P[[q1]] <- list()
                    P[[q1]]$name <- sess[[i]]$cond[[j]]$pmod[[q]]$name
                    P[[q1]]$P    <- sess[[i]]$cond[[j]]$pmod[[q]]$param
                    P[[q1]]$h    <- sess[[i]]$cond[[j]]$pmod[[q]]$poly
                }
            }

            if( length(P) == 0 ) { 
                # P$name <- "none"   ## we want to keep length(P) == 0
                # P$h    <- 0
            }

            U[[j]]$P <- P
        }
        SPM$Sess[[i]]$U <- U

        # user specified regressors
        SPM$Sess[[i]]$C <- list(C=matrix(0,0,0), Cname=character(0))
        if(length(sess$regress) > 0) { 
            C     <- sess$regress[[1]]$val
            Cname <- sess$regress[[1]]$name
            if(length(sess$regress) > 1) {
                for(rr in 2:length(sess$regress)) {
                    C <- cbind(C, sess$regress[[rr]]$val)
                    Cname <- c(Cname, sess$regress[[rr]]$name)
                }
            }
            SPM$Sess[[i]]$C$C    <- C
            SPM$Sess[[i]]$C$name <- Cname
        }

    } # end session

    # factorial design
    if( !is.null(fact) ) {
        SPM$factor <- vector("list", length(fact))
        NC <- length( SPM$Sess[[1]]$U ) # number of conditions
        CheckNC <- 1
        for(i in 1:length(fact)) {
            SPM$factor[[i]]$name   <- fact[[i]]$name
            SPM$factor[[i]]$levels <- fact[[i]]$levels
            CheckNC <- CheckNC * SPM$factor[[i]]$levels
        }
        if(CheckNC != NC) {
            stop("factor do not match conditions")
        }
    } else {
        SPM$factor <- list()
    }

    # Globals
    SPM$xGX$iGXcalc <- global
    SPM$xGX$sGXcalc <- "mean voxel value"
    SPM$xGX$sGMsca  <- "session specific"

    # High pass filter
    SPM$xX$K <- vector("list", length=length(sess))
    for(i in 1:length(sess)) {
        SPM$xX$K[[i]] <- list()
        SPM$xX$K[[i]]$HParam <- sess[[i]]$hpf
    }

    # Autocorrelation
    SPM$xVi$form <- cvi

    # Let SPM configure the design

    # Mask
    if(!is.null(mask)) {
        # TODO!!!
    }

    ##
    ## FILE: SPMROOT/spm_fmri_spm_ui.m
    ##
 
    # get design matrix
    SPM <- spm_fMRI_design(SPM)

    nsess <- length(SPM$nscan)

    # high-pass filter
    K <- vector("list", length=nsess)
    for(i in 1:nsess) {
        K[[i]] <- list(HParam=SPM$xX$K[[i]]$HParam,
                       row=SPM$Sess[[i]]$row,
                       RT=SPM$xY$RT)
    }
    SPM$xX$K <- spm_filter(K)

    # intrinsic autocorrelation (Vi)
    cVi <- tolower(SPM$xVi$form)

    nscan <- SPM$nscan # could be a vector (if multiple sessions)!
    if(is.numeric(cVi)) {
        SPM$xVi$Vi <- spm_Ce(nscan, cVi)
        cVi <- sprintf("AR(%0.1f)", cVi)
    } else if(cVi == "none") {
        SPM$xVi$V <- diag(sum(nscan))
        cVi <- "i.i.d"
    } else { # assume AR(0.2)
        SPM$xVi$Vi <- spm_Ce(nscan, 0.2)
        cVi <- "AR(0.2)"     
    }
    SPM$xVi$form <- cVi


    # at this point, SPM maps the files,
    # checks the orientations, and place the handles
    # in SPM$xY$VY
    #### FIXME: we only use the FIRST SESSION here!!!!! ######
    SPM$xY$VY <- (SPM$SPMR$VY[[1L]][] * SPM$SPMR$VY[[1L]]$scl.slope + 
                                        SPM$SPMR$VY[[1L]]$scl.inter)
    if(verbose) {
        cat("Note: we only use the first session for now\n")
        cat("size VY data: ", round(object.size(SPM$xY$VY)/(1024^2),3), "Mb\n")
    }

    # compute Global variate
    GM <- 100
    
    if(verbose) cat("Calculating globals ... ")
    start.time <- proc.time()[3]
    q <- dim(SPM$xY$VY)[4]
    g <- apply(SPM$xY$VY, MARGIN=4, spm_global)
    #for(i in 1:q) {
    #    g[i] <- spm_global(SPM$xY$VY[,,,i])
    #}
    if(verbose) cat("done. [time = ",(proc.time()[3]-start.time),"]\n",sep="")

    gSF <- 100/g
    if(tolower(SPM$xGX$iGXcalc) == "none") {
        # grand-mean scaling (session/subject-specific)
        # recommended for fMRI, because  the scaling of fMRI data
        # can vary between sessions and could mask regional activations
        # SPM Book p119
        for(i in 1:nsess) {
            gSF[SPM$Sess[[i]]$row] <- GM/mean(g[SPM$Sess[[i]]$row])
        }
    } else {
        # "scaling": different scale per scan
        # to remove the effect of 'gain' for every image (which is known
        # to vary from scan to scan)
    }

  
    # Apply gSF to memory-mapped scalefactors to implement scaling
    q <- dim(SPM$xY$VY)[4]
    for(i in 1:q) {
       SPM$xY$VY[,,,i] <- SPM$xY$VY[,,,i] * gSF[i]
    }

    # place global variates in global structure
    SPM$xGX$rg  <- g
    SPM$xGX$GM  <- GM
    SPM$xGX$gSF <- gSF
 

    # Masking structure automatically set to 80% of mean
    TH <- g * gSF * defaults.mask.thresh

    SPM$xM <- list(T=rep(1, q),
                   TH=TH,
                   I=0,
                   VM=list(),
                   xs=list(Masking="analysis threshold"))

    # Design description - for saving and display
    ntr <- integer(nsess)
    for(i in 1:nsess) {
        ntr[i] <- length(SPM$Sess[[i]]$U)
    }
    SPM$xsDes <- list(
        Basis_functions      = SPM$xBF$name,
        Number_of_sessions   = sprintf('%d',nsess),
        Trials_per_session   = sprintf('%-3d',ntr),
        Interscan_interval   = sprintf('%0.2f {s}',SPM$xY$RT),
        High_pass_Filter     = sprintf('Cutoff: %d {s}',SPM$xX$K[[1]]$HParam),
        Global_calculation   = SPM$xGX$sGXcalc,
        Grand_mean_scaling   = SPM$xGX$sGMsca,
        Global_normalisation = SPM$xGX$iGXcalc)
    

    SPM
}


spmr_fmri_estimate <- function(SPM, method="classical") {

    ##
    ## FILE: SPMROOT/config/spm_run_fmri_est.m
    ##

    start.time <- proc.time()[3]
    verbose <- SPM$SPMR$verbose

    if(tolower(method) == "classical") {
        SPM <- spm_spm(SPM)
        SPM$xVol$M <- SPM$SPMR$VY[[1L]]$qto.xyz
        SPM$xVol$iM <- solve(SPM$xVol$M)

        # Automatically set up contrasts for factor designs
        ## TODO !! ##


    } else if(tolower(method) == "bayesian") {
        stop("bayesian estimation is not implemented yet")
    }

    if(verbose) cat("done. [time = ",(proc.time()[3]-start.time),"]\n",sep="")

    SPM
}

spmr_contrasts <- function(SPM, consess=list(), delete=0) {

    # check/augment consess argument
    stopifnot(is.list(consess) && length(consess) > 0)
    for(i in 1:length(consess)) {
        if(is.null(consess[[i]]$type)) {
            if(is.matrix(consess[[i]]$con) && nrow(consess[[i]]$con) > 1) {
                consess[[i]]$type <- "tcon"
            } else {
                consess[[i]]$type <- "fcon"
            }
        }
        if(is.null(consess[[i]]$sessrep)) {
            # same contrast for several sessions?
            # matlabbatch default is "none" (other options are 
            # "repl", "replna", "sess", and "both")

            # see file SPMROOT/config/spm_cf_con.m
            consess[[i]]$sessrep <- "none"
        }
    } 

    # delete old contrasts entries in SPM?
    if(delete) {
        SPM$xCon <- NULL
        SPM$xCon <- list()
        ## TODO!!
    }

    # get name, STAT, con and sessrep
    for(i in 1:length(consess)) {
        if(consess[[i]]$type == "tcon") {
            name    <- consess[[i]]$name
            STAT    <- "T"
            con     <- consess[[i]]$contrast
            # convert to matrix AND transpose to column vector
            if(!is.matrix(con)) {
                con <- as.matrix(con)
            }
            sessrep <- consess[[i]]$sessrep
        } else if(consess[[i]]$type == "fcon") {
            name    <- consess[[i]]$name
            STAT    <- "F"
            con     <- consess[[i]]$contrast
            # transpose!
            con     <- t(con)
            sessrep <- consess[[i]]$sessrep
        } else {
            stop("consess type can only be 'tcon' or 'fcon'")
        }

        if(sessrep != "none") {
            stop("support for multiple sessions not implemented yet!")
            nsession <- length(SPM$Sess)
            if(sessrep == "repl") {

            } else if(sessrep == "replna") {

            } else if(sessrep == "sess") {

            } else if(sessrep == "both") {

            }
        } else {
            cons = list(con)
            names = list(name)
        }

        # loop over created contrasts
        for(k in 1:length(cons)) {
            # Basic checking of contrast:
            # - right number of elements
            # - spm_sp("isinspp/eachinspp", sX, c)
            c <- spm_conman("ParseCon", cons[[k]], SPM$xX$xKXs, STAT)
           
            # Fill-in the contrast structure
            DxCon <- spm_FcUtil(action="Set", name=names[[k]], STAT=STAT, 
                                set_action="c",  value=c, sX=SPM$xX$xKXs)

            # Append to SPM$xCon
            end <- length(SPM$xCon)
            SPM$xCon[[end + 1]] <- DxCon

            # Compute contrast
            SPM <- spm_contrasts(SPM, length(SPM$xCon))
        }

    } # i in 1:length(consess)

    SPM
}

spmr_write_images <- function(SPM, type="SPM", idx) {

    type <- tolower(type)

    if(type == "beta") {
        if(SPM$SPMR$verbose) cat("writing beta images...")
        # write beta
        for(b in 1:ncol(SPM$xX$X)) {
            filename <- paste(SPM$SPMR$output.dir,
                              sprintf("beta_%04d", b), sep="")
            spmr_write_nifti(V=SPM$Vbeta[[b]], output.file=filename, check=0)
        }
        if(SPM$SPMR$verbose) cat(" done.\n")
    }

    if(type == "resms") {
        if(SPM$SPMR$verbose) cat("writing ResMS image ...")
        filename <- paste(SPM$SPMR$output.dir, "ResMS", sep="")
        spmr_write_nifti(V=(SPM$VResMS), output.file=filename, check=0)
        if(SPM$SPMR$verbose) cat(" done.\n")
    }

    if(type == "con") {
        if(SPM$SPMR$verbose) cat("writing contrast images...")
        for(i in 1:length(SPM$xCon)) {
            filename <- paste(SPM$SPMR$output.dir,
                              sprintf("con_%04d", i), sep="")
            # NOTE: we rescale to get raw ResSS (just like SPM)!
            spmr_write_nifti(V=(SPM$xCon[[i]]$Vcon * SPM$xX$trRV), output.file=filename, check=0)
        }
        if(SPM$SPMR$verbose) cat(" done.\n")
    }

    if(type == "spm") {
        if(SPM$SPMR$verbose) cat("writing SPM{t/F} maps...")
        for(i in 1:length(SPM$xCon)) {
            if(SPM$xCon[[i]]$STAT == "F") {
                filename <- paste(SPM$SPMR$output.dir,
                                  sprintf("spmF_%04d", i), sep="")
            } else if(SPM$xCon[[i]]$STAT == "T") {
                filename <- paste(SPM$SPMR$output.dir,
                                  sprintf("spmT_%04d", i), sep="")
            } else {
                stop("unknown contrast type: ", SPM$xCon[[i]]$type)
            }
            spmr_write_nifti(V=SPM$xCon[[i]]$Vspm, output.file=filename, check=0)
        }
        if(SPM$SPMR$verbose) cat(" done.\n")
    }
}

spmr_results <- function(SPM, 
                         Ic=1L, thresDesc="FWE", u=0.05, k=0, 
                         type="table") {

    # get xSPM
    xSPM <- spm_getSPM(SPM, Ic=Ic, thresDesc=thresDesc, u=u, k=k)

    # Summary list of local maxima for entire volume of interest
    TabDat <- spm_list("List", xSPM)

    TabDat
}
