#include <R.h>
#include <Rinternals.h>

#define TINY 5e-2

SEXP resample_d_1(SEXP Rm, SEXP Rvol,
		  SEXP Rout, SEXP Rgradx, SEXP Rgrady, SEXP Rgradz,
		  SEXP Rx, SEXP Ry, SEXP Rz, 
		  SEXP Rxdim, SEXP Rydim, SEXP Rzdim,
		  SEXP Rbackground, SEXP Rscale, SEXP Roffset) {

    int m = INTEGER(Rm)[0];
    double *vol = REAL(Rvol);
    double *out = REAL(Rout);
    double *gradx = REAL(Rgradx);
    double *grady = REAL(Rgrady);
    double *gradz = REAL(Rgradz);
    double *x = REAL(Rx);
    double *y = REAL(Ry);
    double *z = REAL(Rz);
    int xdim = INTEGER(Rxdim)[0];
    int ydim = INTEGER(Rydim)[0];
    int zdim = INTEGER(Rzdim)[0];
    double background = REAL(Rbackground)[0];
    double *scale = REAL(Rscale);
    double *offset = REAL(Roffset);

    int i;

    for (i=0; i<m; i++)
        {
                double xi,yi,zi;
                xi=x[i]-1.0;
                yi=y[i]-1.0;
                zi=z[i]-1.0;
                if (zi>=-TINY && zi<zdim+TINY-1 &&
                    yi>=-TINY && yi<ydim+TINY-1 &&
                    xi>=-TINY && xi<xdim+TINY-1)
                {
                        double k111,k112,k121,k122,k211,k212,k221,k222;
                        double dx1, dx2, dy1, dy2, dz1, dz2;
                        int off1, off2, offx, offy, offz, xcoord, ycoord, zcoord;

                        xcoord = (int)floor(xi); dx1=xi-xcoord; dx2=1.0-dx1;
                        ycoord = (int)floor(yi); dy1=yi-ycoord; dy2=1.0-dy1;
                        zcoord = (int)floor(zi); dz1=zi-zcoord; dz2=1.0-dz1;

                        xcoord = (xcoord < 0) ? ((offx=0),0) : ((xcoord>=xdim-1) ? ((offx=0),xdim-1) : ((offx=1   ),xcoord));
                        ycoord = (ycoord < 0) ? ((offy=0),0) : ((ycoord>=ydim-1) ? ((offy=0),ydim-1) : ((offy=xdim),ycoord));
                        zcoord = (zcoord < 0) ? ((offz=0),0) : ((zcoord>=zdim-1) ? ((offz=0),zdim-1) : ((offz=1   ),zcoord));

                        /*
                        off1 = xcoord  + xdim*ycoord;
                        off2 = off1+offy;
                        k222 = GET(vol[zcoord     ][off1]); k122 = GET(vol[zcoord     ][off1+offx]);
                        k212 = GET(vol[zcoord     ][off2]); k112 = GET(vol[zcoord     ][off2+offx]);
                        k221 = GET(vol[zcoord+offz][off1]); k121 = GET(vol[zcoord+offz][off1+offx]);
                        k211 = GET(vol[zcoord+offz][off2]); k111 = GET(vol[zcoord+offz][off2+offx]);
                        */
            k222 = vol[(zcoord+0)*(xdim*ydim) + (ycoord+0)*xdim + (xcoord+0)];
            k122 = vol[(zcoord+1)*(xdim*ydim) + (ycoord+0)*xdim + (xcoord+0)];
            k212 = vol[(zcoord+0)*(xdim*ydim) + (ycoord+1)*xdim + (xcoord+0)];
            k112 = vol[(zcoord+1)*(xdim*ydim) + (ycoord+1)*xdim + (xcoord+0)];
            k221 = vol[(zcoord+0)*(xdim*ydim) + (ycoord+0)*xdim + (xcoord+1)];
            k121 = vol[(zcoord+1)*(xdim*ydim) + (ycoord+0)*xdim + (xcoord+1)];
            k211 = vol[(zcoord+0)*(xdim*ydim) + (ycoord+1)*xdim + (xcoord+1)];
            k111 = vol[(zcoord+1)*(xdim*ydim) + (ycoord+1)*xdim + (xcoord+1)];

                        /* WARNING: SPM has different order: z, y and x */
                        gradz[i] = (((k111 - k211)*dy1 + (k121 - k221)*dy2))*dz1 + (((k112 - k212)*dy1 + (k122 - k222)*dy2))*dz2;
                        

                        k111 = (k111*dx1 + k211*dx2);
                        k121 = (k121*dx1 + k221*dx2);
                        k112 = (k112*dx1 + k212*dx2);
                        k122 = (k122*dx1 + k222*dx2);

                        grady[i] = (k111 - k121)*dz1 + (k112 - k122)*dz2;

                        k111 = k111*dy1 + k121*dy2;
                        k112 = k112*dy1 + k122*dy2;

                        gradx[i] = k111 - k112;
                        out[(zcoord)*(xdim*ydim) + (ycoord)*xdim + (xcoord)]   = k111*dz1 + k112*dz2;
                }
                else
                {
                    /*
                        out[(zcoord)*(xdim*ydim) + (ycoord)*xdim + (xcoord)]   = background;
                    */
                        gradx[i] = 0.0;
                        grady[i] = 0.0;
                        gradz[i] = 0.0;
                }
        }

        /* return */
        return R_NilValue;
}

