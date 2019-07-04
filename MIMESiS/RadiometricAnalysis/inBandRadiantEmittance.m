function M = inBandRadiantEmittance(lambda_i,lambda_f,T)
%
% This fucntion computes the IN BAND RADIANT EMITTANCE of a black body at a
% given temperature T under the Lambertian source assumption. This is the
% integral of the radiance with respect a solid angle over the hemisphere
% into which the surface radiates. NOTE: pi and not 2pi due to Lambert's
% cosine law.
%
% INPUT
%   T           temperature [K]
%
% OUTPUT
%   L           radiance [W/m^2]
% 
% -------------------------------------------------------------------------
% Author: Cirelli Renato, Valentina Marchese
% Date: 22/05/2019
% Revision: 1
%
% ChangeLog
% 22/05/2019 - First Version of the file
%
% -------------------------------------------------------------------------
% LICENSED UNDER Creative Commons Attribution-ShareAlike 4.0 International
% License. You should have received a copy of the license along with this
% work. If not, see <http://creativecommons.org/licenses/by-sa/4.0/>.
% -------------------------------------------------------------------------

% Radiant Emittance
M = pi*inBandRadiance(lambda_i,lambda_f,T);

end