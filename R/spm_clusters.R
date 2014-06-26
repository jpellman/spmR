#' @name spm_clusters
#function A = spm_clusters(L)
#' @title SPM: Return the cluster index for a point list
#' @usage [A] = spm_clusters(L)
#' @param L locations [x y x]' {in voxels} ([3 x m] matrix)
#
#' @return cluster index or region number ([1 x m] vector)
#__________________________________________________________________________
#
#' @description spm_clusters characterizes a point list of voxel values defined with
#' their locations (L) in terms of edge, face and vertex connected
#' subsets, returning a list of indices in A, such that the ith location
#' belongs to cluster A(i) (using an 18 connectivity scheme).
#__________________________________________________________________________
# MATLAB version: Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

spm_clusters <- function(L) {

    if(length(L) == 0L) return(numeric(0))

    # Turn location list to binary 3D volume
    #--------------------------------------------------------------------------
    dim    <- apply(L, 2, max)
    vol    <- array(0L, dim=dim)
    vol[L] <- 1L

    # Label each cluster in 3D volume with its own label using an 18 
    # connectivity criterion
    #--------------------------------------------------------------------------
    cci <- spm_bwlabel(vol, 18L)  

    # Map back to list
    #--------------------------------------------------------------------------
    A <- cci[L]
    A
}

