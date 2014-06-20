# spmR: an R package for fMRI data analysis

This is the source code for spmR, an R implementation of Statistical Parametric Mapping, a mode of analysis used in neuroimaging.  

As of ~~8~~ 19 June, 2014, the present maintainer (John Pellman) has obtained the source code from the original author (Yves Rosseel) but ~~has not modified it in any way~~ has modified it only lightly.  By placing the code on GitHub, the present maintainer hopes to make it readily accessible for future development.

The most often used implementation of SPM is written and maintained by the Wellcome Trust at UCL.  The Wellcome Trust's implementation is built upon MATLAB, is reasonably stable, and should be used by most sane people.  The present maintainer's vision for spmR is that it should be used as both an enjoyable diversion and a starting base should neuroscientists ever choose to migrate to the R platform/ecosystem for doing their analyses.

To wit, there are several advantages that an implementation of SPM built upon R might provide:
* Easier access to novel statistical methods which are unlikely to be implemented in MATLAB/Octave (having a user base consisting primarily of statisticians tends to create this advantage).
* R is free (both as in freedom and in beer).
* R readily interfaces with HPC software, such as OpenMPI/MPICH/LAM-MPI.  It also is really friendly with Hadoop.

# Relevant links:
Slides from a presentation by the original author describing spmR: http://www.wias-berlin.de/workshops/neuroimaging2011/slides/rosseel.pdf

The official website of the creators of SPM: http://www.fil.ion.ucl.ac.uk/spm/

# To do:
* Rewrite documentation in comments at the heads of functions to documentation that can be roxygenized.
* Add support for preprocessing.
* Any of the various other TODOs that I come across in the code- these should be added to this list as they are discovered.
* Rewrite documentation for unported spm functions as roxygen comments.
* Go over code and look for potential bugs. 
