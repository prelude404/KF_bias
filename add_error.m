function [z,u,xt,xp] = add_error(x,K,dt)
global imu_bias imu_noise uwb_noise;
u = x(7:9,:);
u = u + imu_bias + imu_noise;

A = [1, 0, 0, dt, 0, 0;
     0, 1, 0, 0, dt, 0;
     0, 0, 1, 0, 0, dt;
     0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 1];
B = [0.5*(dt^2), 0, 0;
     0, 0.5*(dt^2), 0;
     0, 0, 0.5*(dt^2);
     dt, 0,  0;
     0,  dt, 0;
     0,  0, dt];

xp = x(1:6,1);
for j = 1:K-1
    xp(:,j+1) = A * xp(:,j) + B * u(:,j);
end

xt = x(1:6,:);

% z = sqrt(xt(1,:).^2+xt(2,:).^2+xt(3,:).^2);
% the system is motivated by iuput with noise and bias
z = sqrt(x(1,:).^2+x(2,:).^2+x(3,:).^2);
z = z + uwb_noise;

end

