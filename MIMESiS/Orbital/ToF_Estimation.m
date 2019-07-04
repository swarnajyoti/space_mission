%
% ToF estimation for an Hohmann transfer with variable pericenter and fixed
% apocenter @ SC location
%
% The script shows and computes all the possible Hohmann transfers that
% allow the carrier to reach Enceladus' surface.
% The analysis is performed within the two body problem (2BP) framework,
% using the already integrated equation of motion (no ODE solver).
% The spacecraft (SC) motion is assumed to be circular with null
% inclination (it moves on the Ecliptic plane) and a constant altitude of
% h_sc km.
%
% The detaching event is set at 180° of true anomaly, over the SC orbit,
% and it consists of a tangent retrograde impulsive maneuvres. The goal is
% to either fly-by or impact the targeted area (TA).
%
% Moreover, the planet rotation is considered for proper TA tracking.  
%
% The script allows to
%
%   - set the initial longitude of TA @ detaching time
%     NOTE:the analysis
%  Longitude measured counterclockwise, starting from the SC
%     eccentricity vector ( eVers points from left to right in the graph)
%
%   - set the SC orbital altitude
%
%   - set the discretization paramters for the analysis
% -------------------------------------------------------------------------
% Author: Cirelli Renato
% Date: 18/05/2019
% Revision: 4
%
% ChangeLog
% 23/04/2019 - First Version of the file
% 28/04/2019 - Added the spherical cap computation and the nadir-target
%              visibility
% 06/05/2019 - Fixed the LaTex text, added the velocity magnitude plot 
% 18/05/2019 - Added the IGFOV and SWATH profile on a generic descening
%              transfer orbit
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

% PARAMETERS

% Altitude of detachment
h_SC = 220; % [km]

% Target spherical cap surface
coverA = 100; %[m^2]

% Target initial longituted
target_long = deg2rad(180+45); % [rad]

% Visibility angles for the SC2C and C2T
SC2C_minAngle = 0; %[°]
C2T_minAngle = 0; %[°]

% Pericenter of the descending trajectory discretization
periStep = 15;
periI = 1*enceladus.radius;
periF = 0.001*enceladus.radius;

% Time discretization step size
timeStep = 30; %[s]

% -------------------------------------------------------------------------
% Figure Initialization
% -------------------------------------------------------------------------
myFrameHandler = initNewFigure('Moon Centred Intertial Reference Frame');
handler_axes_MCIRF = axes;
axis(1.2*(enceladus.radius+h_SC).*[-1 1 -1 1 -1 1]);
title('Moon Centred Intertial Reference Frame');
daspect([1 1 1]);
view([0,90])
xlabel('$x \;[km]$')
ylabel('$y \;[km]$')
zlabel('$z \;[km]$')
hold on;

initNewFigure('SC - Carrier Distance');
handler_axes_1 = axes;
title('SC - Carrier Distance');
xlabel('$Descent \; Time \;[s]$');
ylabel('$Distance \;[km]$')
hold on;

initNewFigure('Carrier - Target Distace');
handler_axes_2 = axes;
title('Carrier - Target Distace');
xlabel('$Descent \; Time \;[s]$');
ylabel('$Distance \;[km]$')
hold on;

initNewFigure('Carrier - Nadir Distance');
handler_axes_3 = axes;
title('$Carrier - Nadir Distance$');
xlabel('$Descent \; Time \;[s]$');
ylabel('$Distance \;[km]$')
hold on;

initNewFigure('Descending Time vs. Impact velocity');
handler_axes_4 = axes;
title('Descending Time vs. Impact velocity');
xlabel('$Descent \; Time \;[s]$');
ylabel('$Impact \; Velocity \;[m/s]$')
hold on;

initNewFigure('Velocity Variation at Detach');
handler_axes_5 = axes;
title('Velocity Variation at detach');
xlabel('$Descent \; Time \;[s]$');
ylabel('$\Delta V \;[m/s]$')
hold on;

initNewFigure('Velocity Magnitude');
handler_axes_6 = axes;
title('Velocity Magnitude');
xlabel('$Descent \; Time \;[s]$');
ylabel('$\|\underline{V}\| \;[m/s]$')
hold on;

initNewFigure('Velocity Magnitude');
handler_axes_7 = axes;
title('Velocity during the descent');
xlabel('$Local \; Altitude \;[km]$');
ylabel('$\|\underline{V}\| \;[m/s]$')
hold on;

initNewFigure('IGFOV');
handler_axes_8 = axes;
title('IGFOV');
xlabel('$Descent \; Time \;[s]$');
ylabel('$IGFOV \;[m]$')
hold on;

initNewFigure('Swath');
handler_axes_9 = axes;
title('Swath');
xlabel('$Descent \; Time \;[s]$');
ylabel('$Swath \;[m]$')
hold on;
% -------------------------------------------------------------------------
% SC orbit definition
SC_kepVect_I = [enceladus.radius+h_SC,0,0,0,0,pi];
plotOrbitGeo(handler_axes_MCIRF,SC_kepVect_I,[0,2*pi],'b');
plotOrbitGeo(handler_axes_MCIRF,SC_kepVect_I,SC_kepVect_I(6),'bO');
planetPlot(handler_axes_MCIRF,'enceladus');

% SC @ Detachment Event
SC_timeAtDetach =  theta2time(SC_kepVect_I,pi,enceladus.GM);
[SC_xVect_I,SC_vVect_I] = kep2car(SC_kepVect_I,enceladus.GM);

% Spherical Cap angle over the surfae of the moon to achieve coverA m^2 of
% surface
spherCapAngle = acos(1-coverA/(2*pi*enceladus.radius^2));

% Adapted Optic
%myAdaptedOptic [checkLimit,EFL,rad2deg(IFOV),rad2deg(FOV),D_aperture,IGFOV_noadapt,swath_noadapt]
pixelSize = [9e-6,9e-6];
myAdaptedOptic = adaptOptic(220000,pixelSize,30,8000,1.4);

% Pericenter Distance To Test
r_peri_test = linspace(periI,periF,periStep);

for k = 1:length(r_peri_test)
    
    % Compute the Hohmann transfer parameters to define the carrier orbit
    a_car = (norm(SC_xVect_I) + r_peri_test(k))/2;
    e_car = (norm(SC_xVect_I) - r_peri_test(k))/(norm(SC_xVect_I) + r_peri_test(k));
    car_kepVect_I = [a_car,e_car,0,0,0,pi];
    
    % State vector @ detachment
    [car_xVect_I,car_vVect_I] = kep2car(car_kepVect_I,enceladus.GM);
    
    % Compute the impact tru anomaly referred to the carrier orbit
    theta_impact = 2*pi - acos((1/e_car)*(a_car*(1-e_car^2)/enceladus.radius - 1));
    
    % State vector and keplerian vector @ impact event
    car_kepVect_F = [a_car,e_car,0,0,0,theta_impact];
    [car_xVect_F,car_vVect_F] = kep2car(car_kepVect_F,enceladus.GM);
    
    % Compute the time to descent
    car_timeAtDetach = theta2time(car_kepVect_I,pi,enceladus.GM);
    descent_time = theta2time(car_kepVect_I,[pi,theta_impact],enceladus.GM);
    
    % Compute the true anomaly of the SC @ impact event
    SC_thethaAtImpact = time2theta(SC_kepVect_I,SC_timeAtDetach+descent_time,enceladus.GM);
    SC_kepVect_F = [enceladus.radius+h_SC,0,0,0,0,SC_thethaAtImpact];
    [SC_xVect_F,SC_vVect_F] = kep2car(SC_kepVect_F,enceladus.GM);
    
    % Compute the impact location local horizon DCM
    temp = car_kepVect_F(6);
    localHorizDCM = [cos(temp),sin(temp),0;-sin(temp),cos(temp),0;0,0,1];
    handler_impact = plotRefSys3D(handler_axes_MCIRF,car_xVect_F,localHorizDCM,{'n','t','h'},50);

    plotOrbitGeo(handler_axes_MCIRF,car_kepVect_I,[pi,theta_impact],'r--');
    
    % Preallocation fo the cycle
    time_test = 0.01:timeStep:descent_time;
    % Legend
    % 1     SC-Carrier distance
    % 2     Carrier-Target Distance
    % 3     Carrier-Nadir Distance
    % 4     Carrier-Targer Visibility
    % 5     SC-Carrier Visibility
    % 6     Nadir-Target Visibility
    % 7     FOV for full target
    % 8     deltaX at full coverage
    % 9     Velocity Magnitude
    myCollector = zeros(length(time_test),7);
    IGFOV_noadapt = zeros(length(time_test),1);
    swath_noadapt = zeros(length(time_test),1);
    
    for ii = 1:length(time_test)
        
        
        % Current time
        deltaTimeNow = time_test(ii);
        
        % True anomalies at current time on SC and carrier orbit
        SC_thethaNow = time2theta(SC_kepVect_I,SC_timeAtDetach+deltaTimeNow,enceladus.GM);
        car_thethaNow = time2theta(car_kepVect_I,car_timeAtDetach+deltaTimeNow,enceladus.GM);
        
        % Update the keplerian vectors
        SC_kepVectNow = SC_kepVect_I;
        SC_kepVectNow(6) = SC_thethaNow;
        
        car_kepVectNow = car_kepVect_I;
        car_kepVectNow(6) = car_thethaNow;
        
        % Plot the SC and Carrier
        handler_SC = plotOrbitGeo(handler_axes_MCIRF,SC_kepVectNow,SC_kepVectNow(6),'bO');
        handler_car = plotOrbitGeo(handler_axes_MCIRF,car_kepVectNow,car_kepVectNow(6),'rs');
        
        [SC_xVectNow,SC_vVectNow] = kep2car(SC_kepVectNow,enceladus.GM);
        [car_xVectNow,car_vVectNow] = kep2car(car_kepVectNow,enceladus.GM);
        
        % Carrier Velocity
        myCollector(ii,9) = norm(car_vVectNow);
        
        % Compute & Plot SC-Carrier uplink
        myCollector(ii,1) = norm(SC_xVectNow-car_xVectNow);
        temp = [SC_xVectNow;car_xVectNow];
        handler_link1 = plot3(handler_axes_MCIRF,temp(:,1),temp(:,2),temp(:,3),'c:');
 
        % Target Plot
        targetThetaNow = target_long + (enceladus.omega/(60*60*24))*time_test(ii);
        targetDCM = [cos(targetThetaNow),sin(targetThetaNow),0;-sin(targetThetaNow),cos(targetThetaNow),0;0,0,1];
        target_xVectNow = (targetDCM'*[enceladus.radius;0;0])';
        handler_target = plotRefSys3D(handler_axes_MCIRF,target_xVectNow,targetDCM,{'Tz','Te','Tn'},50);
        
        % Compute & Plot Carrier-Target line
        CT_distanceNow = norm(car_xVectNow-target_xVectNow);
        myCollector(ii,2) = CT_distanceNow;
        temp = [car_xVectNow;target_xVectNow];
        handler_link2 = plot3(handler_axes_MCIRF,temp(:,1),temp(:,2),temp(:,3),'b:');
        
        % Carrier-Targer visibility
        myCollector(ii,4) = dot((car_xVectNow -target_xVectNow)./norm((car_xVectNow -target_xVectNow)),targetDCM(1,:));
        
        % Nadir Plot
        rVers_now = car_xVectNow/norm(car_xVectNow);
        nadir_xVectNow = enceladus.radius.*rVers_now;
        nadirDCM = [cos(car_thethaNow),sin(car_thethaNow),0;-sin(car_thethaNow),cos(car_thethaNow),0;0,0,1];
        handler_nadir = plotRefSys3D(handler_axes_MCIRF,nadir_xVectNow,nadirDCM,{'Gx','Gy','Gz'},50);
        
        % Compute & Plot Carrier-Nadir line
        myCollector(ii,3) = norm(car_xVectNow-nadir_xVectNow);
        temp = [car_xVectNow;nadir_xVectNow];
        handler_link3 = plot3(handler_axes_MCIRF,temp(:,1),temp(:,2),temp(:,3),'r:');
        
        % SC-Carrier Visibility
        % Spacecract location in respect the local horizon at narid point
        myCollector(ii,5) = dot((SC_xVectNow -nadir_xVectNow)./norm((SC_xVectNow -nadir_xVectNow)),nadirDCM(1,:));
       
        % Nadir-Target Visibility
        % Check if the narid poitning is inside the target area on hround
        % Angle between nadir and targer
        angleNT = targetThetaNow-car_thethaNow;
        % Check if the nadir poitning is inside the target spherical cap
        if abs(angleNT)< spherCapAngle
            myCollector(ii,6) = 1;
        else
            myCollector(ii,6) = 0;
        end
        
         pause(0.01);
         %drawnow
%         myFrame = getframe(myFrameHandler);
%         myImage = frame2im(myFrame);
%         [imind,cm] = rgb2ind(myImage,256);
%         if ii == 1
%             imwrite(imind,cm,['myAnim_',num2str(k),'.gif'],'gif','Loopcount',inf,'DelayTime',0.1);
%         else
%             imwrite(imind,cm,['myAnim_',num2str(k),'.gif'],'gif','WriteMode','append','DelayTime',0.1);
%         end
        
        % Observed Scene
        % myAdaptedOptic [checkLimit,EFL,(IFOV),(FOV),D_aperture,IGFOV_noadapt,swath_noadapt]
        IGFOV_noadapt(ii) = (pixelSize(1)*1e3*myCollector(ii,2))/myAdaptedOptic(2);
        swath_noadapt(ii) = 2*1e3*myCollector(ii,2)*tan(myAdaptedOptic(4)/2);
        
        delete([handler_SC,handler_car,handler_link1,handler_link2,handler_link3,handler_target,handler_nadir])
        
        
    end
    
    plotOrbitGeo(handler_axes_MCIRF,car_kepVect_I,theta_impact,'r*');
    plotOrbitGeo(handler_axes_MCIRF,SC_kepVect_I,SC_thethaAtImpact,'b*');
    
    delete(handler_impact);
    
    % Plot the SC-Carrier Distance
    temp = myCollector(:,1);
    tempLimit = sin(deg2rad(SC2C_minAngle));
    plot(handler_axes_1,time_test,temp,'r:');
    plot(handler_axes_1,time_test(myCollector(:,5)>=tempLimit),temp(myCollector(:,5)>=tempLimit),'b-');
    
    % Plot the Carrier-Target Distance
    temp = myCollector(:,2);
    tempLimit = sin(deg2rad(C2T_minAngle));
    plot(handler_axes_2,time_test,temp,'r:')
    plot(handler_axes_2,time_test(myCollector(:,4)>=tempLimit),temp(myCollector(:,4)>=tempLimit),'b-')
    
    % Plot the Carrier-Nadir Distance
    temp = myCollector(:,3);
    plot(handler_axes_3,time_test,temp,'r:');
    plot(handler_axes_3,time_test(myCollector(:,6) == 1),temp(myCollector(:,6) == 1),'b-');
    
    % Plot the descent time vs. impact velocity
    plot(handler_axes_4,descent_time,norm(car_vVect_F)*1e3,'bO')
    
    % Plot the deltaV at departure
    plot(handler_axes_5,descent_time,norm(SC_vVect_I-car_vVect_I)*1e3,'bO')
    
    % Plot the deltaV at departure
    plot(handler_axes_6,time_test,myCollector(:,9)*1e3,'b-')
    
    % Plot the deltaV at departure
    plot(handler_axes_7,myCollector(:,3),myCollector(:,9)*1e3,'b-')
    
    % Plot IGFOV prifile (Target Pointing)
    %plot(handler_axes_8,time_test,IGFOV_noadapt,'b-')
    
    temp = IGFOV_noadapt;
    tempLimit = sin(deg2rad(C2T_minAngle));
    plot(handler_axes_8,time_test,temp,'r:')
    plot(handler_axes_8,time_test(myCollector(:,4)>=tempLimit),temp(myCollector(:,4)>=tempLimit),'b-')
    
    
    % Plot the Swath profile )Target Poinitng)
    %plot(handler_axes_9,time_test,swath_noadapt,'b-')
    
    temp = swath_noadapt;
    tempLimit = sin(deg2rad(C2T_minAngle));
    plot(handler_axes_9,time_test,temp,'r:')
    plot(handler_axes_9,time_test(myCollector(:,4)>=tempLimit),temp(myCollector(:,4)>=tempLimit),'b-')
    
    fprintf('\n[%d - %d] \n',k,periStep);
    fprintf('Pericenter position is at %.3f km \n', r_peri_test(k));
    fprintf('deltaV Hohmann = %.3f km/s (%.3f m/s) \n',norm(SC_vVect_I-car_vVect_I),norm(SC_vVect_I-car_vVect_I)*1e3);
    fprintf('Descent Time = %.3f h (%.3f min) \n',descent_time/(60*60),descent_time/(60));
    fprintf('Impact velocity = %.3f km/s (%.3f m/s) \n\n',norm(car_vVect_F),norm(car_vVect_F)*1e3);
    
end