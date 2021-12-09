% clc;
close all;
clearvars -except imu_bias imu_noise uwb_noise K dt sigma_bias sigma_noise sigma_y;
% x0 = [0;0;0;0;0;0]; % start point;
x0 = [50;50;50];
% x0 = [0;0;0;0;0;0];
t = 0:dt:K*dt; % time vars
t(K+1) = [];

% ground truth p_curve,v_curve,u_curve
omega = 10;
% p_curve = [x0(1) + 200*(cos(omega * pi * t/150)-1); % T = 300/omega
%            x0(2) + 0.0000001*(t.^4);
%            x0(3) + 0.0001*(t.^3)];
% v_curve = [-200*(omega*pi/150)*sin(omega * pi * t/150);
%            0.0000004*(t.^3);
%            0.0003*(t.^2)];
% u_curve = [-200*((omega*pi/150)^2)*cos(omega * pi * t/150);
%            0.0000012*(t.^2);
%            0.0006*t];

p_curve = [x0(1)+200*(cos(omega*pi*t/200)-1); %T = 2*pi/(omega*pi/150) = 300/omega;
           x0(2)+200*sin(omega*pi*t/150);
           x0(3)+200*sin(omega*pi*t/100)];  %T = 2*pi/(2*omega*pi/150) = 150/omega;
v_curve = [-200*(omega*pi/200)*sin(omega*pi*t/200);
           200*(omega*pi/150)*cos(omega*pi*t/150);
           200*(omega*pi/100)*cos(omega*pi*t/100)];
u_curve = [-200*((omega*pi/200)^2)*cos(omega * pi * t/200);
           -200*((omega*pi/150)^2)*sin(omega * pi * t/150);
           -200*((omega*pi/100)^2)*sin(omega * pi * t/100)];

% figure(1)
% plot3(p_curve(1,:),p_curve(2,:),p_curve(3,:),'m',x0(1),x0(2),x0(3),'og',0,0,0,'dr');
% xlabel('x');ylabel('y');zlabel('z');
% grid on;

figure(2)
subplot(3,3,1)
plot(t,p_curve(1,:),'r');ylabel('x');
subplot(3,3,2)
plot(t,p_curve(2,:),'r');ylabel('y');
subplot(3,3,3)
plot(t,p_curve(3,:),'r');ylabel('z');
subplot(3,3,4)
plot(t,v_curve(1,:),'r');ylabel('vx');
subplot(3,3,5)
plot(t,v_curve(2,:),'r');ylabel('vy');
subplot(3,3,6)
plot(t,v_curve(3,:),'r');ylabel('vz');
subplot(3,3,7)
plot(t,u_curve(1,:),'r');ylabel('ux');
subplot(3,3,8)
plot(t,u_curve(2,:),'r');ylabel('uy');
subplot(3,3,9)
plot(t,u_curve(3,:),'r');ylabel('uz');

X(1:3,:) = p_curve(:,1:K);
X(4:6,:) = v_curve(:,1:K);
X(7:9,:) = u_curve(:,1:K);

% measurement uwb,imu and x_predict
[uwb, u_real, x_real, x_predict] = add_error(X, K, dt);
y = uwb(:,1:K);
u = u_real(:,1:K);
% [x_KF] = KF_bias(y, u, x_real(:,1), K, dt);
% y_result(y,u,X,K,dt,t);

y_real = sqrt(p_curve(1,1:K).^2 + p_curve(2,1:K).^2 + p_curve(3,1:K).^2);
[x_KF] = KF_bias(y_real, u_curve, u, x_real(:,1), K, dt);
% y_result(y_real, u_curve(:,1:K), X, K, dt,t);

% figure(3)
% plot3(p_curve(1,:),p_curve(2,:),p_curve(3,:),'r');
% hold on
% plot3(x_real(1,:),x_real(2,:),x_real(3,:),'b');
% % plot3(x_KF(1,:),x_KF(2,:),x_KF(3,:),'g');
% legend('real','expected');

figure(4)
[~] = plot_result(t(1:K),x_KF(1:6,:),x_real(1:6,:),'Curve');

figure(5)
% plot(t(1:K),p_curve(1,1:K),'-.k',t(1:K),p_curve(2,1:K),'-.k',t(1:K),p_curve(3,1:K),'-.k','LineWidth',1);
hold on
plot(t(1:K),x_real(1,:),'r',t(1:K),x_real(2,:),'b',t(1:K),x_real(3,:),'m','LineWidth',1);
plot(t(1:K),x_KF(1,:),'--r',t(1:K),x_KF(2,:),'--b',t(1:K),x_KF(3,:),'--m','LineWidth',1);
h1 = legend('x_{gtd}','y_{gtd}','z_{gtd}','x_{esti}','y_{esti}','z_{esti}','Location','northwest','NumColumns',3);
% plot(t,p_curve(1,:),'m',t,p_curve(2,:),'m',t,p_curve(3,:),'m');
% h1 = legend('x_{gtd}','y_{gtd}','z_{gtd}','x_{esti}','y_{esti}','z_{esti}','x_{predict}','y_{predict}','z_{predict}','Location','northwest','NumColumns',3,'FontName','Times New Roman','FontSize',12);
ylabel('Position(m)','FontName','Times New Roman','FontSize',16);
xlabel('Time(s)','FontName','Times New Roman','FontSize',16);
set(h1,'Orientation','horizon','Box','on');
title('Estimation with gtd','FontName','Times New Roman','FontSize',16);

figure(6)
plot(t(1:K),x_real(1,:),'r',t(1:K),x_real(2,:),'b',t(1:K),x_real(3,:),'m','LineWidth',1);
hold on
plot(t(1:K),x_predict(1,1:K),'--r',t(1:K),x_predict(2,1:K),'--b',t(1:K),x_predict(3,1:K),'--m','LineWidth',1);
h1 = legend('x_{gtd}','y_{gtd}','z_{gtd}','x_{predict}','y_{predict}','z_{predict}','Location','northwest','NumColumns',3,'FontName','Times New Roman','FontSize',12);
xlabel('Time','FontName','Times New Roman','FontSize',16);
ylabel('Position','FontName','Times New Roman','FontSize',16);
set(h1,'Orientation','horizon','Box','on');
title('Prediction with gtd','FontName','Times New Roman','FontSize',16);

% figure(7)
% plot3(x_real(1,:),x_real(2,:),x_real(3,:),'b','LineWidth',1);
% hold on
% plot3(x_KF(1,:),x_KF(2,:),x_KF(3,:),'--b','LineWidth',1);
% legend('gtd','esti','FontName','Times New Roman','FontSize',12);
% title('Trajectory','FontName','Times New Roman','FontSize',16);

figure(8)
x_real(7,:) = x_real(1:3,1)' * x_real(4:6,1);
x_real(8,:) = x_real(4:6,1)' * x_real(4:6,1);
error = x_KF - x_real;
error_norm = sqrt(error(1,:).^2 + error(2,:).^2 + error(3,:).^2);
plot(t(1:K),error_norm,'k','linewidth',1);
xlabel('Time','FontName','Times New Roman','FontSize',16);
ylabel('Error','FontName','Times New Roman','FontSize',16);
legend('error','FontName','Times New Roman','FontSize',12);
title('Estimation Error','FontName','Times New Roman','FontSize',16);

figure(9)
set(gcf,'Position',[100,20,600,600]);
subplot(3,1,1)
plot(t,error(1,:),'r',t,error(2,:),'b',t,error(3,:),'k','linewidth',1);
legend('x','y','z','FontName','Times New Roman','FontSize',12);
ylabel('Position Error','FontName','Times New Roman','FontSize',16);
grid on
subplot(3,1,2)
plot(t,error(4,:),'r',t,error(5,:),'b',t,error(6,:),'k','linewidth',1);
legend('v_x','v_y','v_z','FontName','Times New Roman','FontSize',12);
ylabel('Velocity Error','FontName','Times New Roman','FontSize',16);
grid on
subplot(3,1,3)
plot(t,error(7,:),'k','linewidth',1);
hold on
plot(t,error(8,:),'r','linewidth',1);
legend('p_0^Tv_0','v_0^Tv_0','FontName','Times New Roman','FontSize',12);
ylabel('x_7,x_8','FontName','Times New Roman','FontSize',16);
grid on

figure(12)
plot(t,error(7,:),'r','linewidth',1);
hold on
plot(t,error(8,:),'b','linewidth',1);
legend('p_0^Tv_0','v_0^Tv_0','FontName','Times New Roman','FontSize',12);
ylabel('Error of  x_7,x_8','FontName','Times New Roman','FontSize',16);
