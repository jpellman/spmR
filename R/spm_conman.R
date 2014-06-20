#' @name spm_conman
#' @title SPM Contrast Manager
#' @description Contrast manager for spmR.
# FORMAT varargout=spm_conman(varargin)
#       - An embedded callback, multi-function function
#       - For detailed programmers comments,
#         see format specifications in main body of program

# FORMAT [c,I,emsg,imsg,msg] = spm_conman('ParseCon',cstr,X,STAT)
# Contrast weights parser: Catch evaluation errors and invalid contrasts
#' @param cstr string (or cellstr) to evaluate to get contrast(s)
#' @param X design matrix
#' @param STAT 'T' or 'F' (for contrast checking)
#' @param c contrast weights matrix
#' @param I logical validity indicator: indicates which rows of cellstr(cstr) generated valid contrasts which were included in c
# emsg       - cellstr of error messages produced during parsing
# imsg       - cellstr of information messages for valid contrasts
# msg        - cellstr of all messages produced during parsing,
#              one cell per string in cstr

spm_conman <- function(action="ParseCon", cstr, X, STAT="F") {

    action <- tolower(action)

    # only 'ParseCon' for now!
    stopifnot(action == "ParseCon")

    if(action == "ParseCon") {
        p <- spm_SpUtil("size", X, 2)
     
        # FIXME: we assume cstr is a single vector/matrix
        # not a list!!!
        # we assume c is already a proper matrix 
        stopifnot(is.matrix(cstr))
        c <- cstr

        # Evaluate individual lines of contrast matrix input
        # NOTE that c is tranposed!
        for(i in 1:ncol(c)) {
            ## NOTE: we do not allow right-padding (yet)!
            ## contrasts should always have the right size!
            if(length(c[,i]) != p) {
                cat("ncol in design matrix = ", p, "\n")
                cat("contrast line ", i, " = ", c[,i], "\n")
                stop("incorrect number of elements in contrast")
            }

            # check contrast line
            # FIXME: TODO!!
            #spm_SpUtil('isCon', X, t(c))  # 'isCon' is only one
        }
        
        return(c)
    } else

    stop("unknown action argument in spm_conman function")
}

