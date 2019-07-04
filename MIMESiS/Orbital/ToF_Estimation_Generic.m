%
% This script is computing the descending trajectory of an object givend
% the pericenter location (i.e. moon surface) and the former eccentricity.
% Limits of observation are also shown (i.e. green and red circles) as well
% as the descending unit (DU) and target (T) link.
%
% -------------------------------------------------------------------------
% Author: Cirelli Renato
% Date: 05/06/2019
% Revision: 1
%
% ChangeLog
% 05/06/2019 - The script has been changed in order to consider a generic
%              descending orbit with arbitrary high eccentricity
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

% ------------------------------------------------------------------------
% Load Enceladus Data
enceladus = enceladusData();
% ------------------------------------------------------------------------

%% PARAMETERS

% Range Limits
range_max = 220; % [km]
range_min = 0.5; % [km]

% Target spherical cap surface
coverA = 20*20; %[m^2]

% Target initial longituted
target_long = deg2rad(180+90+45); % [rad]

% Visibility angles for the SC2C and C2T
SC2C_minAngle = 0; %[°]
C2T_minAngle = 0; %[°]

% DU descending eccentricity
DU_e = linspace(0.2,0.95,10);

% Time discretization step size
timeStep = 30; %[s]

%% ------------------------------------------------------------------------
% Figure Initialization
% -------------------------------------------------------------------------
myFrameHandler = initNewFigure('Moon Centred Intertial Reference Frame');
handler_axes_MCIRF = axes;
axis(1.2*(enceladus.radius_km+1.1*range_max).*[-1 1 -1 1 -1 1]);
title('Moon Centred Intertial Reference Frame');
daspect([1 1 1]);
view([0,90])
xlabel('$x \;[km]$')
ylabel('$y \;[km]$')
zlabel('$z \;[km]$')
hold on;

%% Orbit Definition

% Planet Plot
planetPlot(handler_axes_MCIRF,'enceladus');

% Apocenter location given the perocenter and the eccentricity
DU_rApo = enceladus.radius_km*(DU_e+1)./(1-DU_e);
DU_a = enceladus.radius_km./(1-DU_e);

% Spherical Cap angle over the surfae of the moon to achieve coverA m^2 of
% surface
spherCapAngle = acos(1-coverA/(2*pi*enceladus.radius_km^2));

for k = 1:length(DU_e)
    
    % Initial DU true anomaly on descending orbit
    DU_theta_I = pi + pi/4;
    
    % Current DU Keplerian Vector pi-2*pi*rand()
    DU_kepVect_I = [DU_a(k),DU_e(k),0,0,0,DU_theta_I];
    
    % State vector @ detachment
    [DU_xVect_I,DU_vVect_I] = kep2car(DU_kepVect_I,enceladus.GM);
    
    % State vector and keplerian vector @ impact event
    DU_kepVect_F = [DU_kepVect_I(1:5),2*pi];
    [DU_xVect_F,DU_vVect_F] = kep2car(DU_kepVect_F,enceladus.GM);
    
    % Compute the time at which the distance is 1.1 the maximum one 
    % (Bisection Method to Find the maximum Range position)
    % Initial Guess
    test_theta = DU_theta_I;
    test_initStep = deg2rad(10);
    test_distance = range_max + 2e3;
   
    while abs(test_distance-range_max)>1
        test_theta_old = test_theta;
        test_theta = test_theta + test_initStep;
        text_kep = [DU_a(k),DU_e(k),0,0,0,test_theta];
        [test_xVect,~] = kep2car(text_kep,enceladus.GM);
        test_target_DCM = [cos(target_long),sin(target_long),0;-sin(target_long),cos(target_long),0;0,0,1];
        test_target_xVect = (test_target_DCM'*[enceladus.radius_km;0;0])';
        test_distanceOld = test_distance;
        test_distance = norm(test_xVect-test_target_xVect);
        
        if test_distance < range_max
           test_initStep = test_initStep/2;
           test_distance = test_distanceOld;
           test_theta = test_theta_old;
        end
    end
    
    temp = test_theta - 10*test_initStep;
    DU_time_t0 = theta2time(DU_kepVect_I,temp,enceladus.GM);
    DU_timeToPeri = theta2time(DU_kepVect_I,[temp,2*pi],enceladus.GM);
    
    % Compute the impact location local horizon DCM
    temp = DU_kepVect_F(6);
    localHorizDCM = [cos(temp),sin(temp),0;-sin(temp),cos(temp),0;0,0,1];
    %handler_impact = plotRefSys3D(handler_axes_MCIRF,DU_xVect_F,localHorizDCM,{'n','t','h'},50);
    
    % Plot the orbit
    handler_orbitDU = plotOrbitGeo(handler_axes_MCIRF,DU_kepVect_I,[DU_theta_I,2*pi],'b--');
    
    % Time vector of the simulation
    dTimeVect = 0:timeStep:DU_timeToPeri;
    
    for tt = 1:length(dTimeVect)

        % Current time
        timeNow = DU_time_t0 + dTimeVect(tt);
        
        % True anomalies at current time on DU orbit
        DU_thethaNow = time2theta(DU_kepVect_I,timeNow,enceladus.GM);
        
        % Update the keplerian vectors
        DU_kepVectNow = DU_kepVect_I;
        DU_kepVectNow(6) = DU_thethaNow;
        
        % DU current state vector
        [DU_xVectNow,DU_vVectNow] = kep2car(DU_kepVectNow,enceladus.GM);
        
        % Plot the DU
        handler_DU_pos = plotOrbitGeo(handler_axes_MCIRF,DU_kepVectNow,DU_kepVectNow(6),'r.');
        handler_DU_pos.MarkerSize = 10;
        
        % Plot the DU minimum and maximum range
%         handler_circle_DU(1) = drawCirle(handler_axes_MCIRF,DU_xVectNow,[0,0,1],220,'green--');
%         handler_circle_DU(2) = drawCirle(handler_axes_MCIRF,DU_xVectNow,[0,0,1],50*0.5,'red--');
        
        
        % Target Plot
        % Update the current angular position referred to the MCIRF
        target_ThetaNow = target_long + (enceladus.omega/(60*60*24))*dTimeVect(tt);
        target_DCM = [cos(target_ThetaNow),sin(target_ThetaNow),0;-sin(target_ThetaNow),cos(target_ThetaNow),0;0,0,1];
        target_xVectNow = (target_DCM'*[enceladus.radius_km;0;0])';
        %handler_target = plotRefSys3D(handler_axes_MCIRF,target_xVectNow,target_DCM,{'Tz','Te','Tn'},50);
        handler_T = plot3(handler_axes_MCIRF,...
                           target_xVectNow(1),...
                           target_xVectNow(2),...
                           target_xVectNow(3),'r.');
        
        % Plot the T minimum and maximum range
        handler_circle_T(1) = drawCirle(handler_axes_MCIRF,target_xVectNow,[0,0,1],220,'green-');
        handler_circle_T(2) = drawCirle(handler_axes_MCIRF,target_xVectNow,[0,0,1],10*0.5,'red-');
        
        % Compute & Plot DU-Target line
        CT_distanceNow = norm(DU_xVectNow-target_xVectNow);
        temp = [DU_xVectNow;target_xVectNow];
        handler_linkDUTar = plot3(handler_axes_MCIRF,temp(:,1),temp(:,2),temp(:,3),'r:');
        
        pause(0.05)

%% Gif Creator ------------------------------------------------------------
%         drawnow
%         myFrame = getframe(myFrameHandler);
%         myImage = frame2im(myFrame);
%         [imind,cm] = rgb2ind(myImage,256);
%         if tt == 1 && k == 1
%             imwrite(imind,cm,'myAnim_test2.gif','gif','Loopcount',inf,'DelayTime',0.1);
%         else
%             imwrite(imind,cm,'myAnim_test2.gif','gif','WriteMode','append','DelayTime',0.1);
%         end
%% ------------------------------------------------------------------------

        delete([handler_DU_pos,handler_linkDUTar,handler_circle_T,handler_T])
        
    end
    
    handler_orbitDU.Color = 'black';
    handler_orbitDU.LineStyle = ':';
    
end