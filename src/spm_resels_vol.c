#include <R.h>
#include <Rinternals.h>


SEXP spm_resels_vol(SEXP Rmask, SEXP Rxdim, SEXP Rydim, SEXP Rzdim,
		    SEXP RFWHM, SEXP RR) {

    int *mask = INTEGER(Rmask);
    int xdim = INTEGER(Rxdim)[0];
    int ydim = INTEGER(Rydim)[0];
    int zdim = INTEGER(Rzdim)[0];
    double *FWHM = REAL(RFWHM);
    double *R = REAL(RR);

    int i,j,k;
    int p=0, ex=0, ey=0, ez=0, fxy=0, fxz=0, fyz=0, c=0;
    double rx, ry, rz;

    for(i=0;i<xdim;i++) {
        for(j=0;j<ydim;j++) {
	    for(k=0;k<zdim;k++) {
                if(mask[(k+0)*(xdim*ydim) + (j+0)*xdim + (i+0)]) {
		    p++;
		    if(mask[(k+0)*(xdim*ydim) + (j+0)*xdim + (i+1)]) ex++;
		    if(mask[(k+0)*(xdim*ydim) + (j+1)*xdim + (i+0)]) ey++;
		    if(mask[(k+1)*(xdim*ydim) + (j+0)*xdim + (i+0)]) ez++;
                    if(mask[(k+0)*(xdim*ydim) + (j+0)*xdim + (i+1)] &&
	               mask[(k+0)*(xdim*ydim) + (j+1)*xdim + (i+0)] &&
                       mask[(k+0)*(xdim*ydim) + (j+1)*xdim + (i+1)]) fxy++;
		    if(mask[(k+0)*(xdim*ydim) + (j+0)*xdim + (i+1)] &&
		       mask[(k+1)*(xdim*ydim) + (j+0)*xdim + (i+0)] &&
                       mask[(k+1)*(xdim*ydim) + (j+0)*xdim + (i+1)]) fxz++;
		    if(mask[(k+0)*(xdim*ydim) + (j+1)*xdim + (i+0)] &&
                       mask[(k+1)*(xdim*ydim) + (j+0)*xdim + (i+0)] &&
                       mask[(k+1)*(xdim*ydim) + (j+1)*xdim + (i+0)]) fyz++;
		    if(mask[(k+0)*(xdim*ydim) + (j+0)*xdim + (i+1)] &&
		       mask[(k+0)*(xdim*ydim) + (j+1)*xdim + (i+0)] &&
		       mask[(k+1)*(xdim*ydim) + (j+0)*xdim + (i+0)] &&
		       mask[(k+0)*(xdim*ydim) + (j+1)*xdim + (i+1)] &&
		       mask[(k+1)*(xdim*ydim) + (j+1)*xdim + (i+0)] &&
		       mask[(k+1)*(xdim*ydim) + (j+0)*xdim + (i+1)] &&
		       mask[(k+1)*(xdim*ydim) + (j+1)*xdim + (i+1)]) c++;
		}
	    }
	}
    }

    rx = 1.0/FWHM[0];
    ry = 1.0/FWHM[1];
    rz = 1.0/FWHM[2];

    /*
    Rprintf("Number of edges:\n");
    Rprintf("ex = %d ey = %d ez = %d\n", ex, ey, ez);
    Rprintf("fxy = %d fxz = %d fyz = %d\n", fxy, fxz, fyz);
    Rprintf("p = %d c = %d\n", p, c);
    */

    R[0] = p-(ex+ey+ez)+(fyz+fxz+fxy)-c;
    R[1] = (ex-fxy-fxz+c)*rx+(ey-fxy-fyz+c)*ry+(ez-fxz-fyz+c)*rz;
    R[2] = (fxy-c)*rx*ry + (fxz-c)*rx*rz + (fyz-c)*ry*rz;
    R[3] = c*rx*ry*rz;

    /* return */
    return R_NilValue;
}
