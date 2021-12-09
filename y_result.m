function [] = y_result(y,u,X,K,dt,t)
%% real
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

%% gtd
y_gtd = zeros(1,K);
y0 = 0.5 * (X(1,1)^2 + X(2,1)^2 + X(3,1)^2);
y = 0.5 * (X(1,:).^2 + X(2,:).^2 + X(3,:).^2);
u = X(7:9,:);
v = zeros(3,K);
d = zeros(3,K);
for i = 1:K
    v(1,i) = dt * trapz(u(1,1:i));
    v(2,i) = dt * trapz(u(2,1:i));
    v(3,i) = dt * trapz(u(3,1:i));
    d(1,i) = dt * trapz(v(1,1:i));
    d(2,i) = dt * trapz(v(2,1:i));
    d(3,i) = dt * trapz(v(3,1:i));
    y_gtd(i) = y(i) - y0 + 0.5 * (d(1,i)^2 + d(2,i)^2 + d(3,i)^2);
end

%% result
figure(10)
subplot(2,1,1)
plot(t,y_,'r',t,y_gtd,'b','LineWidth',1);
legend('y_{gtd}','y_{real}','Location','northwest','FontName','Times New Roman','FontSize',12);
ylabel('Output','FontName','Times New Roman','FontSize',16);
subplot(2,1,2)
error = abs(y_ - y_gtd);
plot(t,error,'k','LineWidth',1);
legend('error','Location','northwest','FontName','Times New Roman','FontSize',12);
ylabel('Output Error','FontName','Times New Roman','FontSize',16);
