# generalized convolution of inputs (U) with basis set (bf)
# FORMAT [X,Xname,Fc] = spm_Volterra(U,bf,V);
# U          -  input structure array
# bf         -  Basis functions
# V          -  [1 or 2] order of Volterra expansion [default = 1]
#
# X          -  Design Matrix
# Xname      -  names of regressors [columns] in X
# Fc(j).i    -  indices pertaining to input i (and interactions)
# Fc(j).name -  names pertaining to input i   (and interactions)
#___________________________________________________________________________
#
# For first order expansions spm_Volterra simply convolves the causes
# (e.g. stick functions) in U.u by the basis functions in bf to create
# a design matrix X.  For second order expansions new entries appear
# in ind, bf and name that correspond to the interaction among the
# original causes. The basis functions for these effects are two dimensional
# and are used to assemble the second order kernel in spm_graph.m.
# Second order effects are computed for only the first column of U.u.

spm_Volterra <- function(U, bf, V=1) {

    if(V == 2) {
        stop("Volterra == 2 is not implemented yet!")
    }

    Fc <- vector("list", length=length(U))
    X  <- matrix(0, nrow=length(U[[1]]$u), ncol=length(U)*ncol(bf))

    for(i in 1:length(U)) {
        ## FIXME!!
        ## we assume (for now) there is only 1 column in U[[i]]$u
        ##
        x <- U[[i]]$u
        idx <- 1:length(x)
        # for each basis function
        for(p in 1:ncol(bf)) {
            z <- convolve(x, rev(bf[,p]), type="open")
            X[,p + (i-1)*ncol(bf)] <- z[idx]
        }
        Fc[[i]]$i <- 1:ncol(bf) + (i-1)*ncol(bf)
        Fc[[i]]$name <- U[[i]]$name
    }

    attr(X, "Fc") <- Fc

    X
}

