 library(spmR)

## 0. get auditory dataset (auditory.nii.gz) 
# wget http://http://www.da.ugent.be/datasets/auditory.nii.gz .

## 1. describe design

# first session/subject
session1 <- list()
session1$scans <- c("auditory.nii.gz")          # full path to 4D fMRI data
session1$nscans <- 84

# first condition
condition1 <- list(name="Condition 1",
                   onset=c(6, 18, 30, 42, 54, 66, 78),
                   duration=6
                  )
session1$cond <- list(condition1)


# 1. specify 1st level
SPM <- spmr_fmri_specify_1stlevel(RT=7, 
                                  units="scans",
                                  #bf.name="hrf (with time and dispersion derivatives)",
                                  bf.name="hrf",
                                  sess=list(session1),
                                  cvi="AR(1)",
                                  output.dir=""
                                 )
# 2. estimate model
SPM <- spmr_fmri_estimate(SPM, method="classical")

# 3 specify and estimate t - contrasts
Tcontrast1 <- list(type="tcon", name="active > rest",
                   contrast=c(1,0) )

Tcontrast2 <- list(type="tcon", name="active < rest",
                   contrast=c(-1,0) )

SPM <- spmr_contrasts(SPM, consess=list(Tcontrast1, Tcontrast2))

# optionally: write out SPM t/F images
# spmr_write_images(SPM, type="SPM")

# 4. inference: compute (corrected) p-values
Table <- spmr_results(SPM, type="table")

# print table
round(Table[,c(1,2,3,5,6,7,9,10,11,12,13,14)], 3)
attr(Table, "footer")


