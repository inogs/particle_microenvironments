function [RDBL wsink] =  calculate_dbl(R0)
%---------------------------------------------
% Input:
% R0 : particle radius (m)
%---------------------------------------------
% Output:
% DBL : thickness of diffusive boundary layer (m)
%       for sinking particle
%       note, here estimates sinking speed too
% wsink: sinking speed of particle (m/s)
%---------------------------------------------
% Calculates diffusive boundary layer thickness from Re and Sch
% using empirical sinking speed and parameters for O2
% Following De Vries et al., 2014, BGS
 eta = 1.17;         % sinking speed exponent 
 cw  = 2.55;         % sinking speed coefficient m^(1-eta)/s
 wsink = cw*(2*R0)^(eta);
% Following Ploug et al. 2002, L&O
 kvw = 1.19e-6;      % Kinematic viscosity at 34 psu and 15C
 Re  = R0*wsink/kvw;
 Sch = 697;          % Schmidt number for O2 at 34 psu and 15C
 Sh = 1 + 0.619 * Re^0.412 * Sch^(1/3);
 RDBL = R0/Sh;

