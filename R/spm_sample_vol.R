# returns voxel values from a memory mapped image - a compiled routine
# FORMAT X = spm_sample_vol(V,x,y,z,hold);
# V      -  is a memory mapped image volume
# x      -  matrix of x coordinates {pixels}
# y      -  matrix of y coordinates {pixels}
# z      -  matrix of z coordinates {pixels}
# hold   -  sets the interpolation method for the resampling.
#           0          Zero-order hold (nearest neighbour).
#           1          First-order hold (trilinear interpolation).
#           2->127     Higher order Lagrange (polynomial) interpolation using
#                      different holds (second-order upwards). 
#          -127 - -1   Different orders of sinc interpolation. 
# X      -  output image
# 
# OR     [X,dX,dY,dZ] = spm_sample_vol(V,x,y,z,hold);
# Similar to above, except that the derivatives in the three orthogonal
# directions are also returned.

# Matlab version: Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

spm_sample_vol <- function(V, x, y, z, hold=1, gradient=TRUE) {
    background <- 0.0

    map <- V
    x <- as.matrix(x); y <- as.matrix(y); z <- as.matrix(z)
    m <- nrow(x)
    n <- ncol(x)

    stopifnot(m == nrow(y), m == nrow(z),
              n == ncol(y), n == ncol(z))
    stopifnot(abs(hold) <= 127)
    
    scale = 1.0; offset = 0.0;
    xdim = dim(V)[1]; ydim = dim(V)[2]; zdim = dim(V)[3]
    if(gradient == FALSE) {
        if(hold == 0) {
            out <- RESAMPLE_0(m=m*n, vol=map, x=x, y=y, z=z, 
                              xdim=xdim, ydim=ydim, zdim=zdim,
                              background=background, 
                              scale=scale, offset=offset)
        } else if(hold == 1) {
            out <- RESAMPLE_1(m=m*n, vol=map, x=x, y=y, z=z, 
                              xdim=xdim, ydim=ydim, zdim=zdim,
                              background=background, 
                              scale=scale, offset=offset)
        } else {
            out <- RESAMPLE_POLY(m=m*n, vol=map, x=x, y=y, z=z,
                                 xdim=xdim, ydim=ydim, zdim=zdim,
                                 hold=hold+1,
                                 background=background,
                                 scale=scale, offset=offset)
        }
        return(out)

    } else if(hold == 0) {
        stop("This wont work for nearest neighbour resampling.\n")

    } else {
        if(hold == 1) {
            #out <- RESAMPLE_D_1(m=m*n, vol=map, x=x, y=y, z=z,
            #                    xdim=xdim, ydim=ydim, zdim=zdim,
            #                    background=background,
            #                    scale=scale, offset=offset)
            #str(out)
            out <- RESAMPLE_D_1.cversion(m=m*n, vol=map, x=x, y=y, z=z,
                                xdim=xdim, ydim=ydim, zdim=zdim,
                                background=background,
                                scale=scale, offset=offset)

            #str(out)
            #stop("for now")
        } else {
            out <- RESAMPLE_D_POLY(m=m*n, vol=map, x=x, y=y, z=z,
                                   xdim=xdim, ydim=ydim, zdim=zdim,
                                   hold=hold+1,
                                   background=background,
                                   scale=scale, offset=offset)
        }
        return(out)
    }
}

RESAMPLE_0 <- function(m, vol, x, y, z, xdim, ydim, zdim,
                       background=NA, scale=1.0, offset=0.0) {
    
    out <- array(NA, dim=dim(vol))
    for(i in 1:m) {
        xcoord <- floor(x[i]+0.5) 
        ycoord <- floor(y[i]+0.5)
        zcoord <- floor(z[i]+0.5)
        if (xcoord>0 && xcoord<=xdim && ycoord>0 &&
            ycoord<=ydim && zcoord>0 && zcoord<=zdim) {
            out[x[i],y[i],z[i]] <- vol[xcoord, ycoord, zcoord]
        } else {
            out[x[i],y[i],z[i]] <- background
        }
    }
    return(out)
}

RESAMPLE_D_1.cversion <- function(m, vol, x, y, z, xdim, ydim, zdim,
                                  background=as.numeric(NA), 
                                  scale=1.0, offset=0.0) {

    # reserve space for return values here
    out <- array(as.numeric(NA), dim=dim(vol))
    gradx <- matrix(as.numeric(NA), m, 1)
    grady <- matrix(as.numeric(NA), m, 1)
    gradz <- matrix(as.numeric(NA), m, 1)

    tmp <- .Call("resample_d_1", as.integer(m), vol, 
                 out, gradx, grady, gradz,
                 as.numeric(x), as.numeric(y), as.numeric(z), 
                 as.integer(xdim), as.integer(ydim), as.integer(zdim),
                 background=background, scale=scale, offset=offset,
                 package="spmR")

    list(out=out, gradx=gradx, grady=grady, gradz=gradz)
}

RESAMPLE_D_1 <- function(m, vol, x, y, z, xdim, ydim, zdim,
                         background=NA, scale=1.0, offset=0.0) {

    TINY <- 5e-2

    out <- array(NA, dim=dim(vol))
    gradx <- matrix(NA, m, 1)
    grady <- matrix(NA, m, 1)
    gradz <- matrix(NA, m, 1)

    for(i in 1:m) {
        if(z[i] >= (1-TINY) && z[i] < (zdim + TINY) &&
           y[i] >= (1-TINY) && y[i] < (ydim + TINY) &&
           x[i] >= (1-TINY) && x[i] < (xdim + TINY)) {
           
            xcoord <- floor(x[i]); dx1 <- x[i]-xcoord; dx2 <- 1.0 - dx1
            ycoord <- floor(y[i]); dy1 <- y[i]-ycoord; dy2 <- 1.0 - dy1
            zcoord <- floor(z[i]); dz1 <- z[i]-zcoord; dz2 <- 1.0 - dz1
            if(xcoord < 1) xcoord <- 1; if(xcoord > xdim) xcoord <- xdim
            if(ycoord < 1) ycoord <- 1; if(ycoord > ydim) ycoord <- ydim
            if(zcoord < 1) zcoord <- 1; if(zcoord > zdim) zcoord <- zdim

            k222 <- vol[xcoord  , ycoord  , zcoord]
            k122 <- vol[xcoord+1, ycoord  , zcoord]
            k212 <- vol[xcoord  , ycoord+1, zcoord]
            k112 <- vol[xcoord+1, ycoord+1, zcoord]
            k221 <- vol[xcoord  , ycoord,   zcoord+1]
            k121 <- vol[xcoord+1, ycoord,   zcoord+1]
            k211 <- vol[xcoord  , ycoord+1, zcoord+1]
            k111 <- vol[xcoord+1, ycoord+1, zcoord+1]
            

            out[x[i],y[i],z[i]] <- ((k111*dx1 + k211*dx2)*dy1 + 
                                    (k121*dx1 + k221*dx2)*dy2)*dz1 + 
                                   ((k112*dx1 + k212*dx2)*dy1 + 
                                    (k122*dx1 + k222*dx2)*dy2)*dz2
            gradx[i,1] <- ((k111     - k211    )*dy1 + 
                           (k121     - k221    )*dy2)*dz1 + 
                          ((k112     - k212    )*dy1 + 
                           (k122     - k222    )*dy2)*dz2
            grady[i,1] <- ((k111*dx1 + k211*dx2) - 
                           (k121*dx1 + k221*dx2))*dz1 + 
                          ((k112*dx1 + k212*dx2) - 
                           (k122*dx1 + k222*dx2))*dz2
            gradz[i,1] <- ((k111*dx1 + k211*dx2)*dy1 + 
                           (k121*dx1 + k221*dx2)*dy2) - 
                          ((k112*dx1 + k212*dx2)*dy1 + 
                           (k122*dx1 + k222*dx2)*dy2)
        
        } else {
            out[x[i],y[i],z[i]] <- background
        }
    }

    out <- list(out=out, gradx=gradx, grady=grady, gradz=gradz)
    return(out)
}
