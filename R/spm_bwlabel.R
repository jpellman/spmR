#function varargout = spm_bwlabel(varargin)
# Connected component labelling in 2D or 3D - a compiled routine
#
# FORMAT: [L,NUM] = spm_bwlabel(bw,conn)
# 
# The calling interface has been modelled on the 
# image processing toolbox routine bwlabel. See 
# the help for that routine for further information.
#
# Input:
# bw:         : (Binary) image to perform labelling on. Can
#               be 2 or 3D.
# conn        : Connectivity criterion. Could be 6(surface)
#               18(edge) or 26(corner). For a 2D bw these
#               correspond to 4, 8 and 8 respectively.
#
# Output:
# L           : Connected component image, i.e. image where
#               each non-zero voxel in bw will have a value
#               corresponding to its label.
# NUM         : Number of connected components in L.
#
#
# The implementation is not recursive (i.e. will not crash for
# large connected components) and is losely based on
# Thurfjell et al. 1992, A new three-dimensional connected
# components labeling algorithm with simultaneous object
# feature extraction capability. CVGIP: Graphical Models 
# and Image Processing 54(4):357-364.
#_______________________________________________________________________
# MATLAB version: Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

spm_bwlabel <- function(bw, conn=18L) {

    #L <- array(0L, dim=dim(bw))
    #NUM <- numeric(1)
    #tmp <- .Call("spm_bwlabel", bw, as.integer(conn), L, NUM)

    L <- spm_bwlabel_joke(bw) # only conn=18??
    
    L
}

spm_bwlabel_joke <- function(vol){

    ## enlarge volume
    dimen <- dim(vol)[1:3]
    dimens <- dimen+2
    vol2 <- array(0,dim=dimens)
    new <- array(0,dim=dimens)
    vol2[2:(dimens[1]-1),2:(dimens[2]-1),2:(dimens[3]-1)] <- vol
    nrset <- 0

    ## find regions

    for(z in 2:(dimens[3]-1)){
      for(y in 2:(dimens[2]-1)){
        for(x in 2:(dimens[1]-1)){

          # if voxel is active, check whether surrounding voxels are active
          if(vol2[x,y,z]==1){
            if(vol2[x-1,y,z]==1||
            vol2[x,y-1,z]==1||
            vol2[x-1,y-1,z]==1||
            vol2[x+1,y,z]==1||
            vol2[x,y+1,z]==1||
            vol2[x+1,y+1,z]==1||
            vol2[x-1,y+1,z]==1||
            vol2[x+1,y-1,z]==1||

            vol2[x,y,z-1]==1||
            vol2[x-1,y,z-1]==1||
            vol2[x,y-1,z-1]==1||
            vol2[x+1,y,z-1]==1||
            vol2[x,y+1,z-1]==1||

            vol2[x,y,z+1]==1||
            vol2[x-1,y,z+1]==1||
            vol2[x,y-1,z+1]==1||
            vol2[x+1,y,z+1]==1||
            vol2[x,y+1,z+1]==1
            ){
              a <- c(new[x-1,y,z],new[x,y-1,z],new[x-1,y-1,z],new[x+1,y,z],new[x,y+1,z],new[x+1,y+1,z],new[x-1,y+1,z],new[x+1,y-1,z],new[x,y,z-1],new[x-1,y,z-1],new[x,y-1,z-1],new[x+1,y,z-1],new[x,y+1,z-1],new[x,y,z+1],new[x-1,y,z+1],new[x,y-1,z+1],new[x+1,y,z+1],new[x,y+1,z+1])



              null <- ifelse(sum(a)>0,1,0)
              b <- unique(a[a!=0])

              # case 1: region is not numbered yet
              if(null==0)
                  {nrset <- nrset+1
                  new[x,y,z] <- nrset}
              # case 2: one surrounding selected voxel is already numbered
              if(length(b)==1)
                  {new[x,y,z] <- b}
              # case 3: multiple surrounding selected voxels are numbered
              if(length(b)>1)
                  {new[x,y,z] <- min(b)
                  for(i in 2:(length(b))){new[new==sort(b)[i]] <- min(b)}
                  }
            }else{nrset <- nrset+1
              new[x,y,z] <- nrset}

          }
        }
      }
    }

    ## relabel (so that all region numbers are sequential

    for(i in 1:max(new)){
      if(length(new[new==i])==0){
        if(length(new[new>i])==0){
          stop
        } else{
          j <- min(new[new>i])
          new[new==j]<-i
        }
      }
    }


    ## shrink volume back to normal dimensions

    new <- new[2:(dim(vol2)[1]-1),2:(dim(vol2)[2]-1),2:(dim(vol2)[3]-1)]

    new
}


