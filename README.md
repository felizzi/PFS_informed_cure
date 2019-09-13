# PFS_informed_cure
There are different files.
The file _Pola_cure_endpoint_direct.R_ performs a fit on the selected endpoints (typically OS, PFSINV - investigator PFS, PFSIRC - IRC assessed PFS) and returns a set of parameters  in the files whose names start with: 

_Parms_mixed_direct_HE_friendly_ -- parameter estimates
_Variance_mixed_direct_HE_friendly_ -- variance on parameter estimates

The file _Pola_informed_cure_endpoints.R_ uses inputs from the estimations obtained with the _direct_ files 
