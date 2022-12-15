%%  This file estimates roll-angle, from horizontal, of a spinning projectile using just the GPS data.

clearvars
close all

[pos_1,pos_2,dt,t,lambda,ang,radius] = gps_data_gen(); % generate GPS data
ang(3,:) = mod(ang(3,:),2*pi);

R_lambda = rotz(rad2deg(lambda));
%radius = 0.075;
g = 9.81;
G = [0;0;g];
%pos_1 = [pos_1(:,1),pos_1];
%gps_vel_1 = gps_vel(pos_1,dt);
%gps_vel_2 = gps_vel(pos_2,dt);
phi_est2(1) = -lambda; % Initialization required!

for count = 1:length(pos_1)-1
    
    if count==1
          phi_est = 0;  % because no velocity estimate yet
%           pos_last_1=pos_1(:,1);
%           pos_last_2=pos_2(:,1);
    else
         gps_pose_1 = [pos_last_1,pos_1(:,count)];
         gps_pose_2 = [pos_last_2,pos_2(:,count)];
         gps_vel_1 = gps_vel(gps_pose_1,dt(count));
         gps_vel_2 = gps_vel(gps_pose_2,dt(count));
        
         vel = 0.5*(gps_vel_1 + gps_vel_2); % averaging vel. of 2 rotating receivers
         theta = atan(-vel(3)/sqrt(vel(1)^2 + vel(2)^2));
         psi = atan(vel(2)/vel(1));
         thetaa(count) = theta;
         psii(count) = psi;

        R_B_NED = [cos(theta)*cos(psi)  -sin(psi) sin(theta)*cos(psi);
                   cos(theta)*sin(psi) cos(psi) sin(theta)*sin(psi);
                   -sin(theta) 0 cos(theta)];     % Body(NR) to NED
        
        rot_mat = R_B_NED*R_lambda*R_B_NED' - eye(3); 

        %c2 = (R_B_NED*R_lambda*R_B_NED' - eye(3))\(pos_1(:,count) - pos_2(:,count));
        %c2 = pinv(rot_mat)*(pos_1(:,count) - pos_2(:,count));
        c2 =  lsqminnorm((rot_mat),(pos_1(:,count) - pos_2(:,count)));
        c1 = R_B_NED*R_lambda*R_B_NED'*c2;

        p = cross(vel,G); % reference horizontal vector
        
    %    orntn_mat = [p,c1,[0;0;1]];
    %    detm = det(orntn_mat);
        
       % phi_est(count) = acos(dot(c2,p)/(norm(c2)*norm(p)));
        n_vec = cross(p,c1)/(norm(c1)*norm(p)); % vector normal to projectile face
        n_vecb = R_B_NED'*n_vec; % in body frame
     
       % num = dot(n_vec,cross(c2,p));
       % den = dot(c2,p);
       %phi_est(count) = atan2(num,den);
        phi_est(count) = acos(dot(p,c1)/(norm(c1)*norm(p)));      %          

        if  n_vecb(3)<0
            phi_est(count) =  2*pi - phi_est(count);%mod(phi_est(count),2*pi);
        else
            phi_est(count) = phi_est(count);
        end
            phi_est2(count) = phi_est(count) - lambda;
            %phi_est(count) = phi_est(count);
    
        %% Simulation of roll estimation
     %{
            viscircles([0 0],radius);
            delete(findall(gcf,'Tag','test'));
            hold on
            plot([0 radius*cos(phi_est(count))], [0 radius*sin(phi_est(count))],'b',[0 radius*cos(phi_est2(count))], [0 radius*sin(phi_est2(count))],'r','Tag','test');
            hold on
            plot([0 radius*cos(ang(3,count)+pi/2)], [0 radius*sin(ang(3,count)+pi/2)],'k','Tag','test');
            legend('GPS-1','GPS-2', 'Dynamics based')
            xlabel(' X (m)')
            ylabel(' Y (m)')
            title('Roll-angle \phi simulation')
            axis equal
            pause(0.5)
    %}
    end
             pos_last_1=pos_1(:,count);
             pos_last_2=pos_2(:,count);
             
                         
end
%% Plots

%plot(gps_vel_1(1,:))
%  plot(t(1:length(thetaa)),rad2deg(psii),'k',t,rad2deg(ang(1,:)),'r')
   figure(2)
%   plot(t(1:length(thetaa)),rad2deg(thetaa),'k',t,rad2deg(ang(2,:)),'r')
subplot(3,1,1)
 %plot(t(1:length(phi_est)),rad2deg(phi_est),'k',t(1:length(phi_est)),rad2deg(ang(3,1:length(phi_est))),'b')
plot(t(1:500),rad2deg(phi_est(1:500)),'k',t(1:500),rad2deg(ang(3,11:510)),'b')
 legend('GPS based', 'Dynamics based')
xlabel('Time (s)')
ylabel('Phi: \phi (\circ)')

subplot(3,1,2)
 plot(t(1:length(phi_est)),rad2deg(thetaa),'k',t(1:length(phi_est)),rad2deg(ang(2,(1:length(phi_est)))),'b')
 legend('GPS based', 'Dynamics based')
xlabel('Time (s)')
ylabel('Theta: \theta (\circ)')

subplot(3,1,3)
 plot(t(1:length(phi_est)),rad2deg(psii),'k',t(1:length(phi_est)),rad2deg(ang(1,(1:length(phi_est)))),'b')
 legend('GPS based', 'Dynamics based')
 xlabel('Time (s)')
ylabel('Psi: \psi (\circ)')
print(figure(2),'angles','-depsc');

figure(3)
plot(t(1:1000),rad2deg(phi_est(1:1000)),'k',t(1:1000),rad2deg(ang(3,11:1010)),'b')
 legend('GPS based', 'Dynamics based')
 xlabel('Time (s)')
ylabel('Phi: \phi (\circ)')
print(figure(3),'roll-angle','-depsc');

% figure(4)
% plot(t(1:end-12),rad2deg(phi_est(1:end-10)-ang(3,1:end-12))-90,'k*')
% ylim([-100 100])
% %plot(t(1:end-2),rad2deg(phi_est(1:end))-rad2deg(ang(3,11:end-2)),'k*')
% %legend('GPS based', 'Dynamics based')
%  xlabel('Time (s)')
% ylabel('Error in \phi (\circ)')
% print(figure(4),'roll-error','-depsc');

