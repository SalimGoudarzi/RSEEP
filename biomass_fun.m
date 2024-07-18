%**************************************************************************
function [DBH] = biomass_fun(age)
%growth curve parameters
beta           = 0.018;
t0             = 6;
dmax           = 119;
%diameter at breast height (DBH) in cm
DBH            = dmax.*(1-exp(-beta*(age-t0)));
DBH(DBH<0)     = 0;

