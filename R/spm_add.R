# add a series of images - a compiled routine
# FORMAT s = spm_add(VI,VO)
# VI    - Vector of mapped volumes (from spm_map or spm_vol).
# VO    - Description of output volume that gets passed to
#         spm_write_plane.m
# flags - Flags can be:
#               'm' - masks the mean to zero or NaN wherever
#                     a zero occurs in the input images.
# s     - Scalefactor for output image.

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

