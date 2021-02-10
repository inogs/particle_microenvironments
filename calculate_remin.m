function Rem0 =  calculate_dbl(R0)
%---------------------------------------------
% Input:
% R0 : particle radius (m)
%---------------------------------------------
% Output:
% Rem0: average remineralization rate in the particle (mmolC/m3/s) 
%       calculated from POC content of particle and specific 
%       carbon remineralization rate
%---------------------------------------------
% Parameters:
 kRemPoc = 1.00/(24*3600);	% 1/s (typical: 1 d-1) specific decay rate of POC in particle
% Particle carbon content from radius R0 and Alldredge et al. (1998) fractal particle parameters
 bPoc = 0.5;                    % 0.5 Fractal exponent for particle mass: M=a*V^b
 aPoc = 1.0;                    % 1.0 ugC/mm^(3b) coefficient for fractal particle mass (note units)
 cnvrs = 10^(-6+9*bPoc+3)/12;   % conversion factor
 Rem0 = cnvrs * (4/3*pi)^(bPoc-1) * aPoc * ...
            kRemPoc * R0^(3*bPoc-3);

