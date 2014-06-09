# computes the number of resels in a volume - a compiled routine
# FORMAT R = spm_resels_vol(V,W)
# V      -  is a memory mapped image volume.
#           Finite and non-zero values are considered to be part of
#           the search volume.
# W      -  smoothness of the component fields {FWHM in voxels}.
# R      - Resel counts, where:
#          R(1) - Euler Characteristic of the volume (number of connected
#                 components - number of holes).
#          R(2) - Resel Diameter (average over all rotations of the
#                 distance between two parallel planes tangent to the
#                 volume in resel space).
#          R(3) - Resel Surface Area (half the surface area of the
#                 volume in resel space).
#          R(4) - Resel Volume (the volume in resel space).
#_______________________________________________________________________
#
# Reference : Worsley KJ et al 1996, Hum Brain Mapp. 4:58-73
#_______________________________________________________________________
# Copyright MATLAB version (C) 2008 Wellcome Trust Centre for Neuroimaging
spm_resels_vol <- function(mask, FWHM) {

    storage.mode(mask) <- "integer"
    R <- numeric(4)
    DIM <- dim(mask); xdim <- DIM[1]; ydim <- DIM[2]; zdim <- DIM[3];
    tmp <- .Call("spm_resels_vol", mask, xdim, ydim, zdim, 
                  FWHM, R, package="spmR")

    R
}

spm_resels_vol_RVERSION <- function(mask, FWHM) {

    xdim <- dim(mask)[1]
    ydim <- dim(mask)[2]
    zdim <- dim(mask)[3]

    Ex <- Ey <- Ez <- Fxy <- Fxz <- Fyz <- C <- 0

    for (i in 1:xdim){
        for (j in 1:ydim){
            for (k in 1:zdim){
                if(mask[i,j,k]==1){
                    Ex <- ifelse(mask[i+1,j,k]==1,Ex+1,Ex)
                    Ey <- ifelse(mask[i,j+1,k]==1,Ey+1,Ey)
                    Ez <- ifelse(mask[i,j,k+1]==1,Ez+1,Ez)

                    Fxy <- ifelse(mask[i+1,j,k]==1 && mask[i,j+1,k]==1 && mask[i+1,j+1,k]==1,Fxy+1,Fxy)
                    Fxz <- ifelse(mask[i+1,j,k]==1 && mask[i,j,k+1]==1 && mask[i+1,j,k+1]==1,Fxz+1,Fxz)
                    Fyz <- ifelse(mask[i,j+1,k]==1 && mask[i,j,k+1]==1 && mask[i,j+1,k+1]==1,Fyz+1,Fyz)

                    C <- ifelse(mask[i+1,j,k]==1 && mask[i,j+1,k]==1 && mask[i+1,j+1,k]==1 && mask[i,j,k+1]==1 && mask[i+1,j,k+1]==1 && mask[i,j+1,k+1]==1 && mask[i+1,j+1,k+1]==1,C+1,C)
                }
            }
        }
    }

    P <- sum(mask)
    rx <- 1/FWHM[1]
    ry <- 1/FWHM[2]
    rz <- 1/FWHM[3]

    R0 <- P-(Ex+Ey+Ez)+(Fyz+Fxz+Fxy)-C
    R1 <- (Ex-Fxy-Fxz+C)*rx+(Ey-Fxy-Fyz+C)*ry+(Ez-Fxz-Fyz+C)*rz
    R2 <- (Fxy-C)*rx*ry + (Fxz-C)*rx*rz + (Fyz-C)*ry*rz
    R3 <- C*rx*ry*rz

    Res <- c(R0,R1,R2,R3)
    Res
}
