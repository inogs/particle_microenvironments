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
NR = 10; 
NL= 6;

Dw = ones(NL,1) * 6.3 * 1e-10; % Check Nadia 

rCd = 2.3*1e-6; % Cd:C stoichiometric ratio of organic matter

Zeu = 100; % Chek Nadia 
Zr  = 200; % length scale over which remineralization slows down

Cr_opt = 0.8;  % Surface ocean specific carbon remineralization rate (used in 1d model) day-1
Cr_min = 0.1;  % Surface ocean specific carbon remineralization rate (used in 1d model) day-1

Cm     = 6.7; 

alpha  = 2.3; 

Ks     = (0.15+0.33)/2; % day-1


R0_vector=linspace(1e-4, 3e-3, NR);
%R0_vector=exp(linspace(log(min),log(max),NR))

POC    = Cm * R0_vector.^alpha;

N0     = 1.0;

beta   = 3.0; 

N      = N0 * R0_vector.^(-beta); 
Depth_profile = [100;200;300;500;750;1000];
O2_profile  = ones(NL,1) * 180.0;
NO3_profile = ones(NL,1) * 8.0;
Cd_profile  =[0.09; 0.14; 0.19; 0.06; 0.22; 0.57];

RD          = zeros(NL,NR);
RS          = zeros(NL,NR);
fOx         = zeros(NL,NR);
fDen        = zeros(NL,NR);
fSul        = zeros(NL,NR);
RDBL_array  = zeros(NL,NR);
Phi_Cd      = zeros(NL,NR);
CdS         = zeros(NL,NR);
Cdtot       = zeros(NL,1);
Cr          = ( Cr_opt - Cr_min) * exp(-(Depth_profile -Zeu)/Zr) + Cr_min; 


for z=1:NL
    for r=1:NR
    disp('=========')
    r
    z
        O2_inf  = O2_profile(z);
        NO3_inf = NO3_profile(z);
        Cd_inf  = Cd_profile(z);
% Set the particle radius (m)
        R0      = R0_vector(r)

% If needed, calculates the sinking speed and DBL thickness
        [RDBL wsink] = calculate_dbl(R0);
  
        RDBL_array(z,r) = RDBL;
        
% If needed, calculates the remineralization rate inside the particle
        Rem0 = calculate_remin(R0); 

% Calculates the denitrifying and sulfidic radiuses, and the particle
% volumetric fractions for oxic, denitrifying and sulfidic metabolisms
        [RD(z,r) RS(z,r) fOx(z,r) fDen(z,r) fSul(z,r)] = calculate_radiuses(R0,RDBL,Rem0,O2_inf,NO3_inf);

% compute Cadmium flux into a particle         

        Phi_Cd(z,r) = 4.0*pi*Dw(z)*R0*R0/RDBL*Cd_inf;

        if RS(z,r) > 0.0

           Jprec = Phi_Cd(z,r) + rCd * Cr(z)* POC(r);

        else

           Jprec = 0.0;   
        
        end

% compute Cds
        CdS(z,r)     = Jprec/Ks;


        Cdp(z,r)     = CdS(z,r) + rCd * POC(r);

    disp('=========')
    end
    Cdtot(z) = dot(Cdp(z,:),N);
end

 

%RD 
%RS 
%fOx 
%fDen
%fSul
%writemtx(RD,"RD.csv")
