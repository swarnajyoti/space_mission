function L = radiance(T)
%
% This fucntion computes the RADIANCE of a black body at a given 
% temperature T. This is the integral of the spectral radiance over the
% entire spectrum [0,+Inf) and the result is the well known
% Stefan-Boltzmann law.
%
% INPUT
%   T           temperature [K]
%
% OUTPUT
%   L           radiance [W/(m^2*sr)]
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

% Radiance
L = ((2*pi^4*Kb^4)/(15*h^3*c^2))*T^4;

end