%
% Radiometric Analysis of MIMESiS Optical System
%
% BAND's CENTRAL WAVELENGTH CONSIDERED
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

%clear
close all
clc

%All the figure are docked in one window
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultFigureWindowStyle','docked');
set(0,'DefaultTextFontSize',12);
set(0,'DefaultAxesFontSize',12);

%% Load Data
% Load Orbital Library
addpath(genpath('myFunctions'))

% Load some properties
enceladus = enceladusData();
sun = sunData();
saturn = saturnData();
QEs = quantumEfficiencyData();

% % Initi Graphs
% %Total power on each channel
% %Photon Flux
% %Channel EQ
% %electrons produced graph

myTitle = 'ECL Radiance';
initNewFigure(['$',myTitle,'$']);
handler_axes_1 = axes;
title(myTitle);
xlabel('$\lambda_c \; [\mu m]$');
ylabel('$Radiance \; [W/(m^2 \; sr \; \mu m)]$')
hold on;

myTitle = 'Total Incoming Flux on a pixel';
initNewFigure(['$',myTitle,'$']);
handler_axes_2 = axes;
title(myTitle);
xlabel('$\lambda_c \; [\mu m]$');
ylabel('$Flux_{pixel} \; [W/m^2]$')
hold on;

myTitle = 'Quantum Efficiency (QE) for VNIR and MIR channels';
initNewFigure(['$',myTitle,'$']);
handler_axes_3 = axes;
title(myTitle);
xlabel('$\lambda_c \; [\mu m]$');
ylabel('$QE(\lambda) \; [e^-/Photon]$')
hold on;

myTitle = 'Electron Production for VNIR and MIR channels';
initNewFigure(myTitle);
handler_axes_4 = axes;
title(myTitle);
xlabel('$\lambda_c \; [\mu m]$');
ylabel('$Flux_{e^-} \; [1/s]$')
hold on;


initNewFigure('SNR');
handler_axes_5 = axes;
title('SNR');
xlabel('$\lambda_c \; [\mu m]$');
ylabel('$SNR \; [dB]$')
hold on;

initNewFigure('NEP');
handler_axes_6 = axes;
title('NEP');
xlabel('$\lambda_c \; [\mu m]$');
ylabel('$NEP \; [W]$')
hold on;

initNewFigure('SEP');
handler_axes_7 = axes;
title('SEP');
xlabel('$\lambda_c \; [\mu m]$');
ylabel('$SEP \; [W]$')
hold on;

%% Sun contributions

% Sun luminosity (or emitted power under BB assumption)
% luminosity_sun = 3.828e26; %[W] <-- Tabulated one ;) 
sun.surface = (4*pi*(sun.radius_km*1e3)^2); %[m^2]
sun.luminosity = pi*spectralRadiance(lambda_c*1e6,sun.T_effective)*sun.surface; %[W/um]

% Sun power flux (or radiant emittance) at enceladus %[w/m^2]
sun.fluxAtMoon = sun.luminosity / (4*pi*enceladus.distance_SunMoon^2);

% Sun power flux (or radiant emittance) at Saturn %[w/m^2]
sun.fluxAtSat = sun.luminosity / (4*pi*saturn.distance_SunSat^2);

%% Saturn constribution
% NOTE: There is a difference between the flux computed with the mean
% surface area and the computed surface area. (the latter is 10 time bigger)

% Saturn luminosity (or emitted power under BB assumption)
saturn.surface = (4*pi*(saturn.radius_km*1e3)^2); %[m^2]
saturn.luminosity = pi*spectralRadiance(lambda_c*1e6,saturn.T_effective)*saturn.surface; %[W/um]

% Saturn flux (or radiant emittance) at enceladus
saturn.fluxAtMoon_emitted = saturn.luminosity / (4*pi*enceladus.distance_SatMoon^2); %[w/(m^2*um)]

%% Enceladus constribution
% Saturn luminosity (or emitted power under BB assumption)
enceladus.flux_emitted = pi*spectralRadiance(lambda_c*1e6,enceladus.T.mean); %[W/m^2]

%% Albedos computation
% 1 - The flux reaching the body BB is divided by the sphere solid angle
%     (i.e. 4pi) in order to obtain the incoming_radiance [w/(m^2*sr)]
% 2 - The reflection happening on BB's surfae is such that it can be
%     considered as an emitting body with known emitted_radiance equal to
%     the albedo*incoming_radiance
% 3 - The emitted_flux is then computed as emitted_radiance*OMEGA where 
%     OMEGA is the solid angle subtended by considered surface SS on the
%     body BB

% Saturn hemispheric solid angle
%OMEGA_saturn = 2*pi;
% Sun Albedo at Saturn %[w/(m^2)]
saturn.albedo_sun = saturn.geometricAlbedo*sun.fluxAtSat;
% Luminosity of sun albedo at saturn %[w]
saturn.albedo_sun_luminosity = saturn.albedo_sun*saturn.surface; 
% Flux (or radiant emittance) at enceladus %[w/m^2]
saturn.fluxAtMoon_reflected = saturn.albedo_sun_luminosity / (4*pi*enceladus.distance_SatMoon^2);

%% How I would do it (Pag 31 of The Book)

%Radiance of the target (Highest case)
radiance_target = (enceladus.geometricAlbedo*(sun.fluxAtMoon +...
                   saturn.fluxAtMoon_emitted +...
                   saturn.fluxAtMoon_reflected) + ...
                   enceladus.flux_emitted)/pi; %[W/(m^2*sr*um)]

% % Normal Case
% radiance_target = (enceladus.geometricAlbedo*(saturn.fluxAtMoon_emitted) + ...
%                    enceladus.flux_emitted)/pi; %[W/(m^2*sr*um)]
               
% % Radiance of the target (Lowest case in eclipses)            
% radiance_target = enceladus.flux_emitted/pi; %[W/(m^2*sr*um)]

% Integration in the selected bandwidth HP: constant radiance within the
% interval
radiance_target_integr = radiance_target.*(2*lambda_delta*1e6); %[W/(sr*m^2)]
               
% Power arriving at one pixel
pixelPower = (pi/4)*(1)*radiance_target_integr*IFOV^2*D_aperture^2; %[W]
% Flux arriving at one pixel
pixelFlux = (pi/4)*(1)*radiance_target_integr/(fNumber^2); %[W/m^2]
% Total power arriving at the sensor
totalPower.visibleChannel = pixelPower(1:3)*n_visibleSpectralPixel;
totalPower.infraredChannel = pixelPower(4:end)*n_infraredSpectralPixel;

%% Photon flux
% Plank's Constant
h = 6.6260693e-34; %[J/s]
% Speed of light
c = 2.99792458e8; %[m/s]
% Total Photons per unit time entering a sensible pixel area
photonFlux = pixelPower./((h*c)./lambda_c); %[photon/s]
%% Quantum Efficiency or Sensitivity Interpolation
% QE is assumed to be constant in the considered spectral channel, with a
% value equal to the central wavelength one
channel_QEs = [QEs.Si_blue_fcn(lambda_c(1)*1e6),...
               QEs.Si_green_fcn(lambda_c(2)*1e6),...
               QEs.Si_red_fcn(lambda_c(3)*1e6),...
               QEs.InSb_fcn(lambda_c(4:end)*1e6)]/100; %[e-/photon]
%% SNR computation
% Amplification factor 
G = 1; %[-]
% Cross talk degradation since light is spread across multiple pixel here
% expressed as 1 (no cross talk) - cross talk percentage
F_mtf.Si = 1-2/100; %[-]
F_mtf.InSb = 1-3/100; %[-]
% Dark current definition NOTE it should be temperature dependant
darkCurrent.Si = 10; %e-/ms @27°C
darkCurrent.InSb = 450; %e-/ms @-70°C ??
% Saturation level definition
e_sat.Si = 48000; %e-
e_sat.InSb = 1e6; %e-
% Read Out noise
e_readOut.Si = 100; %e-
e_readOut.InSb = 100; %e-

% Exposure time
t_exposure.Si = 500e-6; % s
t_exposure.InSb = 500e-6; % s

% Produced electron on each pixel giventhe QE, exposire tume and cross talk
% efficiency and saturation time
F_mtf.Si = F_mtf.Si*ones(1,3);
e_scene.Si = t_exposure.Si*F_mtf.Si.*photonFlux(1:3).*channel_QEs(1:3); %[e-]
t_sat.Si = (e_sat.Si)./(F_mtf.Si.*photonFlux(1:3).*channel_QEs(1:3)); %[s]

F_mtf.InSb = F_mtf.InSb*ones(1,6);
e_scene.InSb = t_exposure.InSb*F_mtf.InSb.*photonFlux(4:end).*channel_QEs(4:end); %[e-]
t_sat.InSb = (e_sat.InSb)./(F_mtf.InSb.*photonFlux(4:end).*channel_QEs(4:end)); %[s]

% Saturation Limit Check
e_scene.Si <= e_sat.Si
e_scene.InSb <= e_sat.InSb

% Dark Current 
e_dark.Si = t_exposure.Si*darkCurrent.Si; %[e-]
e_dark.InSb = t_exposure.InSb*darkCurrent.InSb; %[e-]

% Signal to noise ratio considering
% Ovserved scene, dark curent
SNR.Si =20*log10((G*e_scene.Si)./sqrt((2^0.5)*G*(e_scene.Si+e_dark.Si)));
SNR.InSb =20*log10((G*e_scene.InSb)./sqrt((2^0.5)*G*(e_scene.InSb+e_dark.InSb)));

% NEP (not in standard form)
NEP.Si = ((e_dark.Si*h*c)./(lambda_c(1:3).*channel_QEs(1:3)))*1e3; %[W]
NEP.InSb = ((e_dark.InSb*h*c)./(lambda_c(4:end).*channel_QEs(4:end)))*1e3; %[W]

% NEP equivalent black body temperature
% TBD

% Digital conversion
% Number of bit in the converter
A2D_bit = 12;
% Number of quantized intervals
A2D_quantum = 2^12;
% SEP saturation equivalent power
SEP.Si = ((e_sat.Si*h*c)./(lambda_c(1:3).*channel_QEs(1:3)))*1e3; %[W]
SEP.InSb = ((e_sat.InSb*h*c)./(lambda_c(4:end).*channel_QEs(4:end)))*1e3; %[W]

% SEP equivalent black body temperature
% TBD

% Digital Resolution
digiRes.Si = (SEP.Si-NEP.Si)./A2D_quantum; %[W/bit]
digiRes.InSb = (SEP.InSb-NEP.InSb)./A2D_quantum; %[W/bit]

%%
% Total power on each channel
% Photon Flux
% Channel EQ
% electrons produced graph

% Radiance ECL
temp = [totalPower.visibleChannel,totalPower.infraredChannel];
plot(handler_axes_1,lambda_c,temp,'O-')
ylim(handler_axes_1,[0,1.05*max(temp)]);

% Photon flux per pixel
plot(handler_axes_2,lambda_c,photonFlux,'O-')
ylim(handler_axes_2,[0,1.05*max(photonFlux)]);

% QE Channels
plot(handler_axes_3,lambda_c,channel_QEs,'O-')
ylim(handler_axes_3,[0.2,1.1]);

% Electron Flux on a pixel
temp = [F_mtf.Si.*photonFlux(1:3).*channel_QEs(1:3),F_mtf.InSb.*photonFlux(4:end).*channel_QEs(4:end)];
plot(handler_axes_4,lambda_c,temp','O-');
ylim(handler_axes_4,[0,1.05.*max(temp)]);

% SNR
hTemp_1 = plot(handler_axes_5,lambda_c,[SNR.Si,SNR.InSb]','O-');
hTemp_2 = plot(handler_axes_5,lambda_c,5.*ones(size(lambda_c)),'--','DisplayName','SNR = 5 dB');
ylim(handler_axes_5,[0,50]);
legend(handler_axes_5,hTemp_2)

% NEP
temp = [NEP.Si,NEP.InSb];
plot(handler_axes_6,lambda_c,temp.','O-')
ylim(handler_axes_6,[0,1.05*max(temp)]);
%legend()

% SEP
temp = [SEP.Si,SEP.InSb];
plot(handler_axes_7,lambda_c,temp.','O-')
ylim(handler_axes_7,[0,1.05*max(temp)]);
%legend()
