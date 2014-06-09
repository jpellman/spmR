#function [N,Z,M,A,XYZ] = spm_max(X,L)
# Sizes, maxima and locations of local excursion sets
# FORMAT [N Z M A XYZ] = spm_max(X,L)
# X     - values of 3-D field
# L     - locations [x y z]' {in voxels}
#
# N     - size of region {in voxels)
# Z     - Z values of maxima
# M     - location of maxima {in voxels}
# A     - region number
# XYZ   - cell array of voxel locations
#__________________________________________________________________________
#
# spm_max characterizes a point list of voxel values (X) and their
# locations (L) in terms of edge, face and vertex connected subsets,
# returning a maxima- orientated list:  The value of the ith maximum is
# Z(i) and its location is given by M(:,i). A(i) identifies the ith
# maximum with a region. Region A(i) contains N(i) voxels, whose
# coordinates are in a 3-by-N(i) array in XYZ{i}.
#
# See also: spm_bwlabel.c and spm_clusters.m
#__________________________________________________________________________
# MATLAB version: Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

spm_max <- function(X,L=NULL) {

    stopifnot(is.matrix(L))
    if(nrow(L) == 0L) return(list(N=integer(0), Z=numeric(0), M=integer(0),
                                  A=integer(0), XYZ=L))

    # %-Ensure that L contains exactly integers
    #%--------------------------------------------------------------------------
    stopifnot(is.integer(L))

    #-Turn location list to binary 3D volume
    #--------------------------------------------------------------------------
    dim    <- apply(L, 2, max)
    vol    <- array(0L, dim=dim)
    vol[L] <- 1L

    # Label each cluster with its own label using an 18 connectivity criterion
    # cci = connected components image volume
    #--------------------------------------------------------------------------
    cci <- spm_bwlabel(vol, 18L)
    num <- max(cci)

    # Get size (in no. of voxels) for each connected component
    # ccs = connected component size
    #--------------------------------------------------------------------------
    ccs <- table(cci)[-1]

    # Get indices into L for voxels that are indeed local maxima (using an 18 
    # neighbour criterion)
    #--------------------------------------------------------------------------
    vol[L] <- X
    Lindex <- spm_get_lm(vol, L)

    M <- L[Lindex,]
    Z <- X[Lindex]
    A <- cci[M]
    N <- ccs[A]

    #  Cell array of XYZ locations of voxels in each cluster
    #--------------------------------------------------------------------------
    XYZ = vector("list", length=max(A))
    for(i in 1:max(A)) {
        XYZ[[i]] <- which(cci == i, arr.ind=TRUE)
    }

    list(N=N, Z=Z, M=M, A=A, XYZ=XYZ)
}
