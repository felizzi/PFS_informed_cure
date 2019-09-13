# PFS_informed_cure
Different files perform a different role in estimating cures and informing new extrapolations on different endpoints. 
The file _Pola_cure_endpoint_direct.R_ performs a fit on the selected endpoints (typically OS, PFSINV - investigator PFS, PFSIRC - IRC assessed PFS) and returns a set of parameters  in the files whose names start with: 

_Parms_mixed_direct_HE_friendly_ -- parameter estimates
_Variance_mixed_direct_HE_friendly_ -- variance on parameter estimates

The file _Pola_informed_cure_endpoints.R_ uses inputs from the estimations obtained with the _direct_ files. 

Thes files depend on other files that contain functions that construct background hazard values for the subjects included in the trial. 

One of such files is the _haz_countries.R_ file, which reads a series of mortality tables (extracted from mortality.org) and builds the background hazard rates (annualized and subject specific). 
