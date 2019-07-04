function Lp_lambda = spectralPhotonRadiance(lambda,T)
%
% This fucntion computes the SPECTRAL PHOTON RADIANCE of a black body at a
% given temperature T for the given wavelengths.
%
% INPUT
%   lambda      wavelength vector [um]
%   T           temperature [K]
%
% OUTPUT
%   Lp_lambda   spectral radiance [photon/(s*m^2*sr*um)]
% 
% -------------------------------------------------------------------------
% Author: Cirelli Renato, Valentina Marchese
% Date: 21/05/2019
% Revision: 1
%
% ChangeLog
% 21/05/2019 - First Version of the file
%
% -------------------------------------------------------------------------
% LICENSED UNDER Creative Commons Attribution-ShareAlike 4.0 International
% License. You should have received a copy of the license along with this
% work. If not, see <http://creativecommons.org/licenses/by-sa/4.0/>.
% -------------------------------------------------------------------------

% Plank's Constant
h = 6.6260693e-34; %[Js]
% Speed of light
c = 2.99792458e8; %[m/s]
% Boltzmann's constant
Kb = 1.380658e-23; %[J/K]

% Spectral Photon Radiance
Lp_lambda = (2e18*c./(lambda.^4)).*(exp((1e6*h*c)./(lambda*Kb*T))-1).^(-1);

end