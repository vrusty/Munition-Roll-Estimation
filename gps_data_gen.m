function [r1,r2,dt,t,lambda,ang,radius] = gps_data_gen()

% This function generates position vector of 2 gps receivers placed on the
% periphery of a spinning projectile, spaced apart by angle lambda, from
% body pose generated using dynamics. All in NED frame. 

% Add GPS uncertainty error through "err_bound" variable

close all

data = load('PGK.mat'); % data generated from dynamics' simulation
% u v w x y z p_a p_f q r psi theta phi_a phi_f 
time = load('PGK2.mat');

pos = data.XMat(4:6,:); % NED body-center 'm'
ang = data.XMat(11:13,:);
p_a = data.XMat(7,:); % shell spin rate
 
t = time.T;

dt = diff(t);
%plot(t,p_a);
V = gps_vel(pos,dt); % Body Vel in NED
%V = [V(:,1),V];
lambda = deg2rad(180); % angle between gps receivers
radius = 0.075; % shell radius
%omega_p = 100;

phi_1 = ang(3,1);   % initial GPS-1 angle. Initialization required!
phi_2 = phi_1 - lambda;%2*pi + phi_1 - lambda; % initial GPS-2 angle
psi_last = 0;

err_bound = 0*0.05;                      % 'x' m error in GPS
        

for i = 1:length(dt)
    
    V_n = V(1,i);
    V_e = V(2,i);
    V_d = V(3,i);
    
    theta = atan(-V_d/sqrt(V_n^2 + V_e^2)); % pseudo pitch
   % thetaa(i)=theta;
    psi = atan(V_e/V_n); % pseudo yaw
  %  psii(i) = psi;
    psi_dot = (psi - psi_last)/dt(i);
    
    R_B_NED = [cos(theta)*cos(psi)  -sin(psi) sin(theta)*cos(psi);
               cos(theta)*sin(psi) cos(psi) sin(theta)*sin(psi);
               -sin(theta) 0 cos(theta)];     % Body(No Roll) to NED
    
    rand_1 = (-ones(3,1) + 2*rand(3,1)); % random number vector between -1 & +1
    rand_1 = rand_1/norm(rand_1); 
    eps_1 = err_bound*rand_1; % GPS-1 uncertainty
    
    rand_2 = (-ones(3,1) + 2*rand(3,1));
    rand_2 = rand_2/norm(rand_2); 
    eps_2 = err_bound*rand_2; % GPS-2 uncertainty
    
           
    c1 = [radius*cos(phi_1);radius*sin(phi_1);0]; % GPS-1 pos body
    c2 = [radius*cos(phi_2);radius*sin(phi_2);0]; % GPS-2 pos body

    r1(:,i) = pos(:,i) + R_B_NED*c1 + eps_1 ;    % GPS-1 pos NED
    r2(:,i) = pos(:,i) + R_B_NED*c2 + eps_2;    % GPS-2 pos NED
    
    phi_dot = p_a(i) + psi_dot*sin(theta);  
       
    % Print to see if computed phi matches the data
    %{    
            viscircles([0 0],radius);  
            delete(findall(gcf,'Tag','test'));
            hold on
            plot([0 radius*cos(phi_1)], [0 radius*sin(phi_1)],'b',[0 radius*cos(phi_2)], [0 radius*sin(phi_2)],'r','Tag','test');
%             phi_est(count)
%             pause(1)
            hold on
            plot([0 radius*cos(ang(3,i))], [0 radius*sin(ang(3,i))],'k','Tag','test');
%              ang(3,count)
             pause(1)
     %}        
     phi_1 = phi_1 + phi_dot*dt(i);
    phi_2 = phi_2 + phi_dot*dt(i);
    % phi_1 = mod(phi_1,2*pi);
    % phi_2 = mod(phi_2,2*pi);
    
    psi_last = psi;
%     plot3(r1(1,i)-pos(1,i),r1(2,i)-pos(2,i),r1(3,i)-pos(3,i),'.r')
%     hold on
%     drawnow;
    
%     plot3(t(i),pos(:,i));
%     hold on 
%     plot3(t(i),r1(:,i));
%     hold on 
%     plot3(t(i),r2(:,i));
%     hold on 
%     plot3(t(i),c1,t(i),c2);
%  pause(0.1)
    %     plot(t(i),rad2deg(theta),'.k',t(i),rad2deg(ang(2,i)),'.r')
%     hold on
    %pause(0.1)
%{ 
 viscircles([0.5 0.5],0.5)
 %axis equal
    
 delete(findall(gcf,'Tag','stream'));

 str1 = 'GPS1';
 annotation('textbox',[0.5+c1(1)/2 0.5+c1(2)/2 0.1 0.1],'String',str1,'FitBoxToText','on','Tag','stream');
 str2 = 'GPS2';
 annotation('textbox',[0.5+c2(1)/2 0.5+c2(2)/2 0.1 0.1],'String',str2,'FitBoxToText','on','Tag','stream');

 annotation('arrow',[0.5 (0.75+c1(1))/1.75],[0.5 (0.75+c1(2))/1.75],'Tag','stream');
 annotation('arrow',[0.5 (0.75+c2(1))/1.5],[0.5 (0.75+c2(2))/1.5],'Tag','stream'); 
 drawnow;
 %}
      
end
%r1 = [r1(:,1),r1];
%r2 = [r2(:,1),r2];
%{
%figure(1) % 3d trajectory
% maxx = max(r1(1,:)-pos(1,:));
% maxx
for i =1:length(pos)
grid on
plot(r1(1,:),r2(1,:),pos(1,:),'b');             

plot3(r1(1,i)-pos(1,i),r1(2,i)-pos(2,i),r1(3,i)-pos(3,i),'.r')
dim1 = [abs((r1(1,i)-pos(1,i)))/max(abs(r1(1,:)-pos(1,:))) abs((r1(2,i)-pos(2,i)))/max(abs(r1(2,:)-pos(2,:))) .1 .1]
str1 = 'GPS_1';
annotation('ellipse',dim1);%,'String',str1,'FitBoxToText','on');

pause(0.1)
 hold on
% %grid on
% plot3(r1(1,i),r1(2,i),r1(3,i),'.b')
% %plot(r1(2,:),r2(2,:),pos(2,:),'b'); 
 hold on
% %grid on
plot3(r2(1,i)-pos(1,i),r2(2,i)-pos(2,i),r2(3,i)-pos(3,i),'.r')
dim2 = [abs((r2(1,i)-pos(1,i))/norm((r2(1,i)-pos(1,i)))) abs((r2(2,i)-pos(2,i))/norm((r2(2,i)-pos(2,i)))) .3 .3];
str2 = 'GPS_2';
annotation('textbox',dim2,'String',str2,'FitBoxToText','on');
% plot3(r2(1,i),r2(2,i),r2(3,i),'.k')
%plot(r1(3,:),r2(3,:),pos(3,:),'b'); 
%plot3(xNorth_hyb(1:end-100),yEast_hyb(1:end-100),zUp_hyb(1:end-100),'r');
hold on
end    
legend('Center', 'GPS_1', 'GPS_2')
title('3D trajectory')
xlabel('North (m)')
ylabel('East (m)')
zlabel('Up (m)')
%}
 %plot(t(1:length(thetaa)),rad2deg(thetaa),'k',t,rad2deg(ang(2,:)),'r')
end
