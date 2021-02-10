% Test of microenvironment calculation for an individual particle,
% with default parameters from the literature

% Physical quantities:
% ------------------------------------------------
% Inputs:
% R0 : particle radius (m)
% RDBL: diffusive boundary layer (DBL) thickness (m) 
%        this can be estimated by "calculate_dbl.m"
% wsink: particle sinking speed (m/s)
%        (used for DBL calculation only)
%        this can be estimated by "calculate_dbl.m"
% Rem0 : POC remineralization rate (mmol/m3/s)
%        this can be estimated by "calculate_remin.m"
% O2_inf: seawater O2 concentration (mmol/m3)
% NO3_inf: seawater NO3 concentration (mmol/m3)
% ------------------------------------------------
% Outputs:
% RD : denitrification radius (remineralization switch from O2 to NO3) (m)
% RS : sulfidic radius (remineralization switch from NO3 to SO4) (m)
% fOx : fraction of particle volume where remineralization uses O2
% fDen : fraction of particle volume where remineralization uses NO3
% fSul : fraction of particle volume where remineralization uses SO4
% ------------------------------------------------

% Set the particle radius (m)
 R0 = 1/1000;

% If needed, calculates the sinking speed and DBL thickness
 [RDBL wsink] = calculate_dbl(R0);

% If needed, calculates the remineralization rate inside the particle
 Rem0 = calculate_remin(R0); 

% Calculates the denitrifying and sulfidic radiuses, and the particle
% volumetric fractions for oxic, denitrifying and sulfidic metabolisms
 O2_inf = 50;
 NO3_inf = 10;
 [RD RS fOx fDen fSul] = calculate_radiuses(R0,RDBL,Rem0,O2_inf,NO3_inf)

