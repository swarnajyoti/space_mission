% MIMESiS THERMAL SYSTEM DATA

% First attempt of one single node simulation. Hypotesis of SSC: 
% Qinternal-Qoutput+Qinput=0
%--------------------------------------------------------------------------
% DONE:
% Check power dissipation from cpu and electric unit
% Check the Saturn magnetosphere heat flux generation
% Check the REF
% Check alpha and eps
% Check IR radiation from Enceladus and Saturn
% Check the flux exchange due to the atmosphere: DeltaT = 30K
% Last version: 21/05/2019

clear all; clc; close all
time = 7200;
set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaulttextInterpreter','latex');
set(0,'defaultTextFontSize',12);
set(0,'DefaultLegendFontSize',5);
set(0,'defaultAxesFontSize',12);


% Temperature distribution plot
%Temperatures needed

Tds = 4;        % T deep space [K]

Y = [238,343    %PCU 290.5
     273,343    %OBDH Memory 308
     233,358    %OBDH Processor 295.5
     157,252];  %Average range T of detectors 204.5
 T_sc = mean(Y,[1 4]);
 
 Tsc_h = T_sc(1,2); %Tsc_h  T high internal spacecraft [K]
 Tsc_c = T_sc(1,1); %Tsc_c  T cold internal spacecraft [K]

 
 Ym = mean(Y); T0 = mean(Ym);
 %      253,353    % Silicon detector(visible) 303K average
%      60,150];   % Indium antimonide detector(infrared) 105K average
figure ('Name','Operational T')
title('Operational Temperature ranges');
b=bar(Y,'BaseValue',T0);bl = b.BaseLine;c1 = bl.Color;
bl.Color = 'blue';
ylabel('T [K]');somenames={'PCDU';'OBDH Memory';'OBDH Processor';...
       'Detector'};%;'Indium antimonide detector'};
set(gca,'xticklabel',somenames) 
ylim([50 700])
grid on;


%% Enceladus' hot case

 % Carrier data

A_front = 1 * 0.1^2;    %Area in front of the planet 1dm^2 [m^2]
A_side = 2 * A_front;   %Lateral area: 2U-cubesat considered [m^2]

Gsc =10;                % Thermal conductiviy [W/m^2K]
    
REF = 1;                % Radiative Effective Factor (changes with altitude)     

Qdiss = 3;              % Total heat dissipated by OBDH and EPS [W]


alpha_f =1;             %Absorbivity of frontal face
eps_f = 1;              %Emissivity of frontal face
alpha_l = 0.87 ;        %7075 alluminium alloy lateral absorbivity 
eps_l =0.81;            %7075 alluminium alloy lateral emissivity

alpha_b =0.87 ;          
eps_b = 0.81;   

% Heat fluxes

Na = 6.022e23;          % Avogadro number [1/mol]
mass_H20 =18.01528*1e-3;% H20 atomic mass [kg/mol]
n_H20 = 12.2e11;        % Numerical value of number density [1/m^3]
rho = n_H20*mass_H20/Na;% Enceladus maximum density [kg/m^3] 
v = 90;                 % Maximum velocity reached [m/s]
R = A_side/(2*0.1);     % Mean radius
Qball = sqrt(rho/R)*v^3;% Peak heating due to Enceladus atmosphere

sigma = 5.67e-8;        % Stefan-B constant
a = 0.99;               % Enceladus Albedo
asat = 0.47;            % Saturn's albedo

L = 3.9e26;             % Luminosity of the Sun [W]
d = 1.496e11 +1.272e12; % 1AU + Enceladus dinstance to Earth [m]         
WS = L / (4*pi*d^2);    % Solar flux density reaching Enceladus [W/m^2]

r = (252.1+200)*1e3;    % Maximum altitude  [m]
Wir = 4.7e9 /(4*pi*r^2);% Heat flux radiation from Enceladus [W/m^2]

r = (237948)*1e3;       % Enceladus dinstance from Saturn [m]
Wsat =5.4;              % Heat flux radiation from Saturn [W/m^2]

% Hot case: Qi+QIR+Qa+QS+QIRsat+Qasat-Qds = 0
Qi = Qdiss;
QS = WS *A_side * REF;
Qa = WS *A_front *alpha_b*a*REF;
Qir  = Wir*A_front *eps_b*REF;

Qirsat = Wsat * A_side *eps_b*REF;
Qasat = WS *A_side *alpha_b*asat*REF;

%% Earth's hot case

au = 1.496e11;
WS_e = L /(4*pi*au^2);  % Solar flux density from 1 AU [W/m^2]
Wir_e = 224;            % Heat flux radiation from Earth [W/m^2]
a_e = 0.3;              % Earth's albedo

QS_e = WS_e *A_front * REF;
Qa_e = WS_e *A_front *alpha_b*a_e*REF;
Qir_e = Wir_e *A_front *eps_b*REF;

model = 'OneNHotCase.slx';
load_system(model)
sim(model)


%% Plot hottest situation

figure ('Name','Hot case 1SN')
        
        Tsc_h_vect = Tsc_h*ones(time,1);
        Tsc_c_vect = Tsc_c*ones(time,1);

        legend();hold on
        title('Transient of Temperature in 1 node hot steady state case');
        plot(Tsc_h_vect ,'LineWidth',2,'DisplayName','Highest T admissible')
        plot(Tsc_c_vect ,'LineWidth',2,'DisplayName','Lowest T admissible')
        plot(Thot_earth ,'LineWidth',2,'DisplayName','T transient for Earth hot case')
        xlabel('time (s)');ylabel('T [K]');
        grid on; grid minor

        
%% Cold case: Qi+QIR-Qd=0

Qi=0.2*Qdiss;
QS=0;
Qa=0;
Qasat=0;
Qirsat=0;

QS_e = 0;
Qa_e = 0;

model = 'OneNColdCase.slx';
load_system(model)
sim(model)


figure ('Name','Cold case 1SN')
        
        legend();hold on
        title('Transient of Temperature in Enceladus cold steady state case')
        plot(Tsc_h_vect ,'LineWidth',2,'DisplayName','Highest T admissible')
        plot(Tsc_c_vect ,'LineWidth',2,'DisplayName','Lowest T admissible')
        plot(Tcold_enc ,'LineWidth',2,'DisplayName','T transient for Enceladus cold case')
        xlabel('time (s)');ylabel('T [K]');
        grid on; grid minor
        
%% OneNode choosing TCS procedure

% alpha_f =1;% 0.14;      % Pag 128 multilayer composition material
% eps_f = 1;%0.6; 
% 
% alpha_l = 0.87 ;        % Pag 87 tesi
% eps_l =0.81;% 0.88;  
% 
% alpha_b =0.87 ;          % Pag 128 multilayer composition material
% eps_b = 0.81;%0.059;