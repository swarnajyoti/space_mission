%
% Performance of an optical system adapted at different altitude of
% observation under Pushbroom Scannin Technique.
%
% This script compute the performance of an imaging system capable to
% obsere a scene from a certain known distance given
% 
% - The Distance from the scene
% - The sensor pixel dimension
% - The ground resolution IGFOV [m]
% - The covered SWATH [m]
% - The optical f-# the system must achieve
% - The spectral band of interest
%
% The produced information are related to the equivalent optic package
% optical properties
%
% - EFL
% - FOV
% - IFOV
% - Aperture diameter
%
% Moreover, the diffraction limit of the resulting optical system is
% considered and different F-# are shown in order to adress the feasibility
% of the measurement.
% 
% -------------------------------------------------------------------------
% Author: Cirelli Renato, Valentina Marchese
% Date: 14/05/2019
% Revision: 1
%
% ChangeLog
% 14/05/2019 - First Version of the file
%
% -------------------------------------------------------------------------
% LICENSED UNDER Creative Commons Attribution-ShareAlike 4.0 International
% License. You should have received a copy of the license along with this
% work. If not, see <http://creativecommons.org/licenses/by-sa/4.0/>.
% -------------------------------------------------------------------------

clear
close all
clc

%All the figure are docked in one window
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultFigureWindowStyle','docked');
set(0,'DefaultTextFontSize',12);
set(0,'DefaultAxesFontSize',12);

% Load Orbital Library
addpath(genpath('myFunctions'))

initNewFigure('Diffraction Limit Check');
handler_axes_1 = axes;
title('Diffraction Limit Check');
xlabel('$\lambda \; [\mu m]$');
ylabel('$Aperture \; Diameter \; [m]$')
hold on;

%% Altitude Adaptation
% Carrier's Aperture Altitude
h = 220000; %[m]
pixel_size_space = 9e-6; %[m]
pixel_size_spectral = 9e-6; %[m]
IGFOV =30; %[m]
swath = 8000;% %[m]
fNumber = 1.4;

n_visibleSpectralPixel = 1;
n_infraredSpectralPixel = 6;

% Spectral Channel Central Frequency 
lambda_c = [[485.5,532.5,685].*1e-3,1.04,1.25,1.5,1.65,2,4.6].*1e-6; %[m]

% Spectral Channel Width
lambda_delta = [[9.5,37.5,65].*1e-3,0.01,0.02,0.1,0.01,0.3,0.005].*1e-6; %[m]

%% SENSOR SIZE
% Number of pixel along the spatial dimension
n_pixel_space = swath/IGFOV;
% Sensor spatial dimension
sensor_size_space = pixel_size_space*n_pixel_space;
% Sensor spectral dimension
% Pushbrum Scanning
sensor_size_spectral.visible = pixel_size_spectral*n_visibleSpectralPixel;
sensor_size_spectral.infrared = pixel_size_spectral*n_infraredSpectralPixel;

%% Scene Parameter
% Field of View given the swath
FOV = 2*atan2(swath/2,h);

% Required focal length
EFL = (pixel_size_space*h)/IGFOV;

% Instantaneous FOV
IFOV = 2*atan2(IGFOV/2,h);

% Aperure size
D_aperture = EFL/fNumber;

% Depth of Focus
% Admissible variation of the detector relative to the nominal focal plane
% position
DOF = 2*(pixel_size_spectral)*fNumber;

%% Behavior at different altitude
h_lowest = 500; %[m]
IGFOV_noadapt = (pixel_size_space*h_lowest)/EFL;
swath_noadapt = 2*h_lowest*tan(FOV/2);
FOV_adapted = 2*atan2(swath/2,h_lowest);

%% Frequency of interest

% Diffraction limit assessed by Rayleigh criterio
D_limit = @(lambda) (1.22*EFL/IFOV).*lambda;

% Diffraction limit check
plot(handler_axes_1,1e6.*lambda_c,D_limit(lambda_c),'b-O','DisplayName','Limit')
hold on

plot(handler_axes_1,1e6.*lambda_c',D_aperture.*ones(length(lambda_c),1),'--','DisplayName',['f-',num2str(fNumber)]);

%% Output Text

fprintf('Optical System Adapted at %.0f km of altitude \n\n',h);
fprintf('EFL: %.3f mm \n', EFL*1e3);
fprintf('DOF: %.3f mm \n', DOF*1e3);
fprintf('FOV: %.3e rad \n', FOV);
fprintf('IFOV: %.3e rad \n', IFOV);
fprintf('Aperture Diameter: %.3f mm \n\n', D_aperture*1e3);

fprintf('Optical System Performance at %.0f km of altitude \n\n',h_lowest);
fprintf('IGFOV: %.3e m \n', IGFOV_noadapt);
fprintf('SWATH: %.3e m \n', swath_noadapt);
fprintf('FOV for constant SWATH: %.3e rad\n',FOV_adapted);
