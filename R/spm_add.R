#' @name spm_add
#' @title SPM: add a series of images
#' @usage s = spm_add(VI,VO)
#' @param VI Vector of mapped volumes (from spm_map or spm_vol).
#' @param VO Description of output volume that gets passed to
#'         spm_write_plane.m
#' flags - Flags can be:
#'               'm' - masks the mean to zero or NaN wherever
#'                     a zero occurs in the input images.
#' @return Scalefactor for output image.

# NOTE: compiled routine in SPM


spm_add <- function(VI) {

    VO <- VI[[1]]

    if( length(VI) > 1) {
        for(i in 2:length(VI)) {
            VO <- VO + VI[[2]]
        }
    }

    VO
}

