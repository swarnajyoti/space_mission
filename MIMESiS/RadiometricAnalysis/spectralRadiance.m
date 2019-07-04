function L_lambda = spectralRadiance(lambda,T)
%
% This fucntion computes the SPECTRAL RADIANCE of a black body at a given
% temperature T for the given wavelengths.
%
% INPUT
%   lambda      wavelength vector [um]
%   T           temperature [K]
%
% OUTPUT
%   L_lambda    spectral radiance [W/(m^2*sr*um)]
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

% Spectral Radiance
L_lambda = (2e24*(h*c^2)./(lambda.^5)).*(exp((1e6*h*c)./(lambda*Kb*T))-1).^(-1);

end