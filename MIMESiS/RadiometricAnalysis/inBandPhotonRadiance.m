function Lp = inBandPhotonRadiance(lambda_i,lambda_f,T)
%
% This fucntion computes the IN BAND PHOTONRADIANCE of a black body at a
% given temperature T. This is the integral of the spectral photon radiance
% over the entire spectrum [lambda_i,lambda_f] by computing the one-side
% integration and then use it to reconstruct the fnite one. This method has
% been introduced by Widger and Woodall and it is converging to at least 10
% digits.
%
% INPUT
%   lambda_i    initial wavelength [um]
%   lambda_f    final wavelength [um]
%   T           temperature [K]
%
% OUTPUT
%   Lp          radiance [photon/(s*m^2*sr)]
%
% -------------------------------------------------------------------------
% REFERENCES
% Widger, W. K. and Woodall, M. P., Integration of the Planck
% blackbody radiation function, Bulletin of the Am. Meteorological Society,
% 57, 10, 1217-1219, Oct. 1976
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

% In-band Rdiance computation
Lp = oneSideIntegration(h,c,Kb,lambda_f,T) - oneSideIntegration(h,c,Kb,lambda_i,T);

end

function Lp_oneside = oneSideIntegration(h,c,Kb,lambda,T)
% One-side integration of the spectral photon radiance from lambda_i to Inf
% in [photon/(s*m^2*sr)]. Note that the integration is done by an iterative
% process and the used variable is the Wavenumber associated to the given
% wavelength [um]

    % Frequency definition [Hz]
    nu = (1e6*c)/lambda;
    % Wavenumber definition [1/cm]
    sigma = nu/(100*c);

    % Initialization
    Lp_oneside = 0;
    x = (100*h*c*sigma)/(Kb*T);
    % Summation
    for k = 1:min(2+20/x,512)
        Lp_oneside = Lp_oneside + ((x^2)/k + (2*x)/k^2 + 2/k^3)*exp(-k*x);
    end
    % Pre-Summation factor
    Lp_oneside = ((2*Kb^3*T^3)/(h^3*c^2))*Lp_oneside;
    
end