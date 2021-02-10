 function [RD RS fOx fDen fSul] = radius_continuous_function(R0,RDBL,Rem0,O2_INF,NO3_INF);
%---------------------------------------------
% Input:
% R0 : particle radius (m)
% Rem0 : POC remineralization rate (mmol/m3/s)
% RDBL : external effective diffusive boundary layer (m)
% O2_INF : O2 concentration in seawater (mmol/m3)
% NO3_INF : NO3 concentration in seawater (mmol/m3)
%---------------------------------------------
% Output:
% RD : denitrification radius (remineralization switch from O2 to NO3) (m)
% RS : sulfidic radius (remineralization switch from NO3 to SO4) (m)
% fOx : fraction of particle volume where remineralization uses O2
% fDen : fraction of particle volume where remineralization uses NO3
% fSul : fraction of particle volume where remineralization uses SO4
%---------------------------------------------
% Parameters (could be set as additional input):
% Diffusivities in seawater and aggregates
 Dw = 1.7e-9;               % molecular diffusion of O2, NO3 in seawater (15C 34 psu) (m2/s)
 Dp = 1.7e-9;               % molecular diffusion of O2, NO3 inside aggregates (m2/s)
% Stoichiometry of remineralization (ignoring NH4+ oxidation within particles)
 rOC = 1/1;                 % O:C for respiration of organic matter
 rNCDen = 4/5;              % NO3:C for denitrification
% Reduction rates
 RemO2   = Rem0 .* rOC;
 RemNO3  = Rem0 .* rNCDen;
% Concentrations for the transition to linear limitation of oxidation:
 KO2R  = 1.0;		% O2 (mmol/m3)
 KNO3R = 1.0;		% NO3 (mmol/m3)
%---------------------------------------------
% Calculates radiuses
%---------------------------------------------
% Some useful quantities
 DDE = Dw/Dp/RDBL;
 R02 = R0^2;
 R03 = R0^3;
%---------------------------------------------
% Radius of oxic/denitrifying transition
%---------------------------------------------
 rho = (Dp*KO2R/RemO2)^0.5;
 RDP = RemO2/Dp;
 % useful functions
 A = @(r) (KO2R*r*(1/3*r^2/rho^2+1-r/rho*coth(r/rho))); 
 B = @(r) (KO2R*r/rho*(coth(r/rho)-r/rho/2));
 C = @(r) (RDP/3*R0 - A(r)/R0^2 - DDE*(O2_INF-RDP/6*R0^2-A(r)/R0-B(r)));
 % Checked against numerical model
 try
    RD = fzero(@(r) C(r),R0/2);
 catch
    RD = 0;
 end
 % For small O2_INF when KO2<O2(R0) the solution should give RD>R0
 % In this case sets RD=R0 to have a valid solution
 O20 = RDP/6*R02+A(RD)/R0+B(RD);
 if O20<=KO2R & RD>R0
    disp(['Setting RD=R0']);
    RD = R0;
 end
 % Set physical limits - this would happen if the algorithm were to fail
 RD(RD<0)  = 0;
 RD(RD>R0) = 0;
%---------------------------------------------
% Radius of oxic/denitrifying transition
%---------------------------------------------
 rho = (Dp*KNO3R/RemNO3)^0.5;
 RDP = RemNO3/Dp;
 RD3 = RD^3;
 % useful functions
 A = @(r) (KNO3R*r*(1/3*r^2/rho^2+1-r/rho*coth(r/rho))); 
 B = @(r) (KNO3R*r/rho*(coth(r/rho)-r/rho/2));
 C = @(r) (RDP/3*RD3/R02 - A(r)/R02 - DDE*(NO3_INF+RDP/3*RD3 * ...
           (1/R0-3/2*1/RD)-A(r)/R0-B(r)));
 try
    RS = fzero(@(r) C(r),RD/2);
 catch
    RS = 0;
 end
 % For small NO3_INF when KNO3<NO3(RD) the solution should give RS>RD
 % In this case sets RD=R0 to have a valid solution
 NO3S = RDP/6*RD^2+A(RS)/RD+B(RS);
 if NO3S<=KNO3R & RS>RD
    disp(['Setting RS=RD']);
    RS = RD;
 end
 % Set physical limits - this would happen if the algorithm were to fail
 RS(RS<0)  = 0;
 RS(RS>RD) = 0;
%---------------------------------------------
% Calculate remineralization fractions from radii:
% assumes remineralization fractions are proportional to the
% volume of each remineralization shell
 if ~isnan(RD)
    fOx  = ((R03-RD^3)/R03);
    if ~isnan(RS)
       fDen = ((RD^3-RS^3)/R03);
       fSul = (RS^3/R03);
    else
       fDen = 1-fOx;
       fSul = 0;
    end
 else
    fOx  = 1;
    fDen = 0;
    fSul = 0;
 end

