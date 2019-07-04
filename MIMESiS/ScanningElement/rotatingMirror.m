clc
close all

% Mirror rotation angles in respect the optical axis

% evalVect = linspace(deg2rad(-10),deg2rad(10),100);
% n_eval = length(evalVect);

% collector_total_x = zeros(n_eval,n_eval);
% collector_total_y = zeros(n_eval,n_eval);
% collector_SFOV = zeros(n_eval,n_eval);

% for ii = 1:n_eval
%     for jj = 1:n_eval

        % Minimum Mirror Angle (Optical Axis Before Reflection as Ref)
        alpha_min = pi/4 - deg2rad(15.5);%evalVect(ii);
        % Maximum Mirror Angle
        alpha_max = pi/4 + deg2rad(22);%evalVect(jj);
        
%         if alpha_min>alpha_max
%             collector_total_x(ii,jj) = nan;
%             collector_total_y(ii,jj) = nan;
%             collector_SFOV(ii,jj) = nan;
%             continue
%         end

        % Minimum Reflected Optical Axis Angle (Nadir Pointing as Ref)
        theta_min = 2*alpha_min - pi/2;
        % Maximum Reflected Optical Axis Angle 
        theta_max = 2*alpha_max - pi/2;

        % SwipeFOV angle (i.e. Pointing Capability of the Mirror)
        if theta_min < 0
            SFOV = abs(theta_max) + abs(theta_min);
        else
            SFOV = abs(theta_max) - abs(theta_min);
        end
%         collector_SFOV(ii,jj) = SFOV;
        
        % Mirror Length to reflect at aplha_min
        L_mirr = D_aperture/sin(alpha_min);
        
        % Space
        mirr_ingrombro_y = 0.5*L_mirr*cos(pi/2-alpha_max);
        mirr_ingrombro_x = 0.5*L_mirr*sin(pi/2-alpha_max);
        
        % Additional space needed to have the mirror @vertical position
        % starting from the first lens edje
        s3 = (L_mirr-D_aperture)/2;
        
        % RHS Obstacle distance from the first lens edje
        s2 = 28/1000;
        
        % Mirror Center of Roation - RHS Obstacle
        s1_star = (0.5*D_aperture)*tan(pi/2 - alpha_max);
        s1 = tan(theta_max)*(D_aperture+s2)+s1_star;
        
%         collector_total_x(ii,jj) = s1 + mirr_ingrombro_x;
%         collector_total_y(ii,jj) = max([(2*mirr_ingrombro_y),(D_aperture+2*s2)]);
        total_dim_x = s1 + 1.1*mirr_ingrombro_x;
        total_dim_y = max([(2*mirr_ingrombro_y),(D_aperture+2*s2)]);
%         
%     end
% end
% 
% temp_alpha_min = 45 + rad2deg(evalVect);
% temp_alpha_max = 45 + rad2deg(evalVect);

% figure(1)
% surf(temp_alpha_min,temp_alpha_max,collector_total_x)
% xlabel('Alpha Min');
% ylabel('Alpha Max');
% zlabel('Ingombro x');
% figure(2)
% surf(temp_alpha_min,temp_alpha_max,collector_total_y)
% xlabel('Alpha Min');
% ylabel('Alpha Max');
% zlabel('Ingombro y');
% figure(3)
% surf(temp_alpha_min,temp_alpha_max,rad2deg(collector_SFOV))
% xlabel('Alpha Min');
% ylabel('Alpha Max');
% zlabel('SFOV');


fprintf('SwipeFOV: %.3f ° \n', rad2deg(SFOV));
fprintf('D_aperture: %.3e cm \n', D_aperture*1e2);
fprintf('L mirror: %.3e cm \n', L_mirr*1e2);
fprintf('Mirror overall size x at alphamax: %.3e cm \n', mirr_ingrombro_x*1e2);
fprintf('Mirror overall size along y at alphamax: %.3e cm \n', mirr_ingrombro_y*1e2);
fprintf('Overall sapce required in front of the first lens (along x): %.3e cm \n', total_dim_x*1e2);
fprintf('Overall sapce required in front of the first lens (along y): %.3e cm \n', total_dim_y*1e2);