#' @name spm_defaults
#' @title SPM: Sets the defaults which are used by SPM
#
#' @usage spm_defaults()
#_______________________________________________________________________
#
#' @description This file is intended to be customised for the site.
#' Individual users can make copies which can be stored in their own
#' matlab subdirectories. If ~/matlab is ahead of the SPM directory
#' in the MATLABPATH, then the users own personal defaults are used.


# Stats defaults
#=======================================================================
defaults.stats.maxmem   = 2^26;
defaults.stats.maxres   = 64;
defaults.stats.fmri.ufp = 0.001;
defaults.stats.pet.ufp  = 0.05;
defaults.stats.eeg.ufp  = 1;
defaults.stats.topoFDR  = 1;

# fMRI design defaults
#=======================================================================
defaults.stats.fmri.fmri_t    = 16;
defaults.stats.fmri.fmri_t0   = 1;
defaults.stats.fmri.cond.tmod = 0;
defaults.stats.fmri.hpf       = 128;
defaults.stats.fmri.hrf.derivs = c(0,0);
defaults.stats.fmri.volt      = 1;
defaults.stats.fmri.global    = "None";
defaults.stats.fmri.cvi       = "AR(1)";

# Factorial design defaults
#=======================================================================
defaults.stats.fact.dept     = 0;
defaults.stats.fact.variance = 1;
defaults.stats.fact.t2.gmsca = 0;
defaults.stats.fact.ancova   = 0;
defaults.stats.fact.mcov.iCC = 1;
defaults.stats.fact.iCFI     = 1;
defaults.stats.fact.iCC      = 1;
defaults.stats.fact.athresh  = 100;
defaults.stats.fact.rthresh  = .8;
defaults.stats.fact.imask    = 1;
defaults.stats.fact.gmsca    = 50;
defaults.stats.fact.glonorm  = 1;
defaults.stats.fact.mreg_int = 1;

# Model estimation defaults
#=======================================================================
defaults.stats.est.signal = "UGL";
defaults.stats.est.ARP    = 3;

# Contrast manager batch defaults
#=======================================================================
defaults.stats.con.delete = 0;

# Results report batch defaults
#=======================================================================
defaults.stats.results.threshtype = "FWE"; # Threshold type
defaults.stats.results.thresh     = 0.05;  # Threshold value
defaults.stats.results.extent     = 0;     # Spatial extent
defaults.stats.results.maskthresh = 0.05;  # (Uncorrected) Threshold for masking
defaults.stats.results.print      = TRUE;  # Print report to file

# Mask defaults
#=======================================================================
defaults.mask.thresh    = 0.8;


