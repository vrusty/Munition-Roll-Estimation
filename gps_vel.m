function vel = gps_vel(pos,dt)
        
% Calculate GPS velocity in NED using GPS pose in NED

        dv_N = diff(pos(1,:));
        %dv_N = [dv_N(1),dv_N];
       % size(dv_N)
        %size(dt)
        v_N = dv_N./dt;
        v_E = diff(pos(2,:))./dt;
        v_D = diff(pos(3,:))./dt;
        vel = [v_N;v_E;v_D];
        
end