#include <R.h>
#include <Rinternals.h>


SEXP spm_global(SEXP Rvol, SEXP Rnvoxels, SEXP RGX) {

    double *vol = REAL(Rvol);
    double *GX = REAL(RGX);
    int nvoxels = INTEGER(Rnvoxels)[0];

    double tmp, mean, cutoff;
    int i, m;

    /* phase 1: compute mean over all voxels V */
    m = 0; tmp = 0.0;
    for(i=0;i<nvoxels;i++) {
        if(R_FINITE(vol[i])) {
            tmp += vol[i];
	    m++;
	}
    }
    mean = tmp/(double)m;

    /* determine cutoff: mean/8 */
    cutoff = mean/8.0;

    /* phase 2: compute mean over all voxels that exceed the cutoff */
    m = 0; tmp = 0.0;
    for(i=0;i<nvoxels;i++) {
	if(R_FINITE(vol[i]) && vol[i] > cutoff) {
            tmp += vol[i];
	    m++;
	}
    }

    (*GX) = tmp/(double)m;

    /* return */
    return R_NilValue;
}
