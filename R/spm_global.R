# returns the global mean for a memory mapped volume image - a compiled routine
# FORMAT GX = spm_global(V)
# V   - memory mapped volume
# GX  - mean global activity
#_______________________________________________________________________
#
# spm_global returns the mean counts integrated over all the  
# slices from the volume
#
# The mean is estimated after discounting voxels outside the object
# using a criteria of greater than > (global mean)/8

# in SPM this a compiled routine: spm_global.c

spm_global <- function(V) {

    #cutoff <- mean(V)/8
    #idx <- which(V < cutoff)
    #V2 <- V[-idx]
    #GX <- mean(V2)
    #cat("GX R version = ", GX, " -- ")

    #V <- as.numeric(V); nvox = length(V)
    GX <- 0.0
    tmp <- .Call("spm_global", V, length(V), GX)
    #cat("GX C version = ", GX, "\n")

    GX 
}

