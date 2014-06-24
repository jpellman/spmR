#' @name spm_get_lm
#function varargout = spm_get_lm(varargin)
#' @title SPM: Identification of local maxima in 3(or 2)D volume - a compiled routine
#
#' @usage INDEX = spm_get_lm(VOL,LIST)
#
#' @description Routine that identifies which voxels in a list of coordinates
#' that are local maxima, and returns a list of indicies into
#' the coordinate list for those maxima.
#
# Input:
#' @param VOL 3(or 2)D volume of statistics (e.g. t or F)
#' @param LIST 3xn (or 2xn) list of voxel coordinates of 
#                tentative local maxima.
#
# Output:
#' @return INDEX        : Index into LIST such that LIST(:,INDEX)
#'                returns those coordinates that are truly
#'                local maxima.
#_______________________________________________________________
# MATLAB version: Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

spm_get_lm <- function(VOL, LIST) {

    INDEX <- spm_get_lm_joke(VOL, LIST)
    idx <- integer(0)
    for(i in 1:nrow(LIST)) { 
        if(!is.na(INDEX[LIST[i,1], LIST[i,2], LIST[i,3]])) 
            idx <- c(idx, i) 
    }

    idx
}

spm_get_lm_joke <- function(vol,cci){

    # enlargen dataset by 1 

    vol[is.na(vol)] <- 0
    di <- dim(vol)
    dimx <- di[1]+2
    dimy <- di[2]+2
    dimz <- di[3]+2
    larger <- array(0,dim=c(dimx,dimy,dimz))
    larger[2:(dimx-1),2:(dimy-1),2:(dimz-1)] <- vol

    # find local maxima in larger 

    lmax <- array(NA,dim=dim(larger)[1:3])
    for(z in 2:(dim(larger)[3]-1)){
      for(y in 2:(dim(larger)[2]-1)){
    for(x in 2:(dim(larger)[1]-1)){
      if(larger[x,y,z]!=0){
        a <- c(larger[x,y,z],larger[x-1,y,z],larger[x+1,y,z],larger[x,y+1,z],larger[x,y-1,z],larger[x,y,z-1],larger[x,y,z+1],larger[x-1,y-1,z],larger[x+1,y-1,z],larger[x-1,y+1,z],larger[x+1,y+1,z],larger[x-1,y,z-1],larger[x-1,y,z+1],larger[x+1,y,z-1],larger[x+1,y,z+1],larger[x,y-1,z-1],larger[x,y+1,z-1],larger[x,y-1,z+1],larger[x,y+1,z+1])
        maxa <- max(a)
        lmax[x,y,z] <- ifelse(larger[x,y,z]==maxa,larger[x,y,z],NA)
      }
    }
      }
    }

    lmax[2:(dim(larger)[1]-1),2:(dim(larger)[2]-1),2:(dim(larger)[3]-1)]
}


