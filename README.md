SFU Guest lecture covering very basic principles of ray sound transmission and ray tracing. 
Uses modified bellhop model in R initially prublisehd 

Raytrace_TL Traces the ray of a sound through a varying soundspeed profile for a fixed amount of time. 
Also plots the provided sound speed profile and all traces generated. 

Original R code is from Taiki Sakai (NWFSC) and available on github https://github.com/TaikiSan21/PAMmisc/blob/137f9e885eb36421c119fdd28c5c6a94f9ccab16/R/raytrace.R#L4
Retrieved March 30th, 2024.  His code is based on written by Val Schmidt from the University of New Hampshire Val Schmidt (2021). raytrace
https://www.mathworks.com/matlabcentral/fileexchange/26253-raytrace), 
MATLAB Central File Exchange. Retrieved June 29, 2021.

This code has been modified to include a very basic TL value of 18 log(r) where r is the distance travelled by the ray. This is an incorrect assumption
but is sufficient for demonstration purposes.

Similarly compute_avg_TL_grid uses the raytract_tl function to generage averaged TL grid from an initial source locaiton. 
