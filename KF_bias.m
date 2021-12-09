function [x_esti] = KF_bias(y,u,imu, x0,K,dt)
%% y_,x0;
y_ = zeros(1,K);
y0 = 0.5 * (y(1).^2);
y = 0.5 * (y.*y);
% integral: u to v to d
v = zeros(3,K);
d = zeros(3,K);
for i = 1:K
    v(1,i) = dt * trapz(u(1,1:i));
    v(2,i) = dt * trapz(u(2,1:i));
    v(3,i) = dt * trapz(u(3,1:i));
    d(1,i) = dt * trapz(v(1,1:i));
    d(2,i) = dt * trapz(v(2,1:i));
    d(3,i) = dt * trapz(v(3,1:i));
    % is trapz ok? or x=x0+vt+0.5at2?
    y_(i) = y(i) - y0 + 0.5 * (d(1,i)^2 + d(2,i)^2 + d(3,i)^2);
end

%% x_k+1 = A * x_k + B * u
A = [1,0,0,dt,0,0,0,0;
     0,1,0,0,dt,0,0,0;
     0,0,1,0,0,dt,0,0;
     0,0,0, 1,0,0,0,0;
     0,0,0, 0,1,0,0,0;
     0,0,0, 0,0,1,0,0;
     0,0,0, 0,0,0,1,0;
     0,0,0, 0,0,0,0,1];
 
B = [0.5*dt^2,0,0;
     0,0.5*dt^2,0;
     0,0,0.5*dt^2;
     dt,   0,   0;
     0,   dt,   0;
     0,    0,  dt;
     0,    0,   0;
     0,    0,   0];

%% initialize
% x0 = x0(1:6,1)+[50;20;20;0;0;0]; % introduce some initial error
x0 = x0(1:6,1) + [0;0;10;0;0;0];
x_esti = zeros(8,K);
x_esti(:,1) = [x0(1:3,1); x0(4:6,1); x0(1:3,1)'*x0(4:6,1); x0(4:6,1)'*x0(4:6,1)];
% x_esti(:,1) = x_esti(:,1) + [6.8711;-295.4539;-31.8324;-0.5322;-2.2126;2.9324;-101.2608; -347.3815];
x_pre = zeros(8,K);
x_pre(:,1) = x_esti(:,1);
L = zeros(8,8);

%% parameter
% P = diag([1e-2,1e-2,1e-2,1e-3,1e-3,1e-3,1e-1,1e-1]);
% Bv = diag([1e-5,1e-5,1e-5,1e-7,1e-7,1e-7,1e-1,1e-4]);
% V = diag([1e-3,1e-3,1e-3,1e-8,1e-8,1e-8,1e-2,1e-4]);
% R = 1e-4;

% limited error of x7,x8
% Bv = diag([1e-5,1e-5,1e-5,1e-7,1e-7,1e-7,1e0,1e-1]);
% V = diag([1e-3,1e-3,1e-3,1e-8,1e-8,1e-8,1e-2,1e-3]);

% output without error
P = diag([1e1,1e1,1e1,1e-1,1e-1,1e-1,1e2,1e2]);
Bv = diag([1e-5,1e-5,1e-5,1e-6,1e-6,1e-6,1e0,1e-4]);
V = diag([1e-2,1e-2,1e-2,1e-4,1e-4,1e-4,1e1,1e-3]);
R = 1e2;


P = diag([1e1,1e1,1e1,1e-1,1e-1,1e-1,1e4,1e4]);
Bv = diag([1e-5,1e-5,1e-5,1e-6,1e-6,1e-6,1e-8,1e-8])*0;
V = diag([1e-2,1e-2,1e-2,1e-4,1e-4,1e-4,1e-8,1e-8]);


for i = 2:K
    x_pre(:,i) = A * x_esti(:,i-1) + B * imu(:,i-1);
    P = A * P * A' + Bv + V + A * L * Bv + Bv * L' * A';
    H = [d(1,i),d(2,i),d(3,i),0,0,0,((i-1)*dt),0.5*(((i-1)*dt)^2)];
    W = H * P * H' + R;
    K = P * H' / W;
    x_esti(:,i) = x_pre(:,i) + K * (y_(:,i) - H * x_pre(:,i));
    P = P - K * W * K';
    L = (eye(8) - K * H) * (A * L + eye(8));
end
