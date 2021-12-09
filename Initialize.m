% run initialize.m first, then curve.m
clc;
clear;
dt = 1/25;
K = 300/dt;
sigma_bias = diag([0.0001,0.0001,0.0001])*0;
sigma_noise = diag([0.001,0.001,0.001]);
sigma_y = 0.001;

% sigma_bias = diag([0,0,0]);
% sigma_noise = diag([0,0,0]);
% sigma_y = 0;

global imu_bias imu_noise uwb_noise;

imu_bias = sqrt(sigma_bias)*randn(3,1);
imu_noise = sqrt(sigma_noise)*randn(3,K);
uwb_noise = sqrt(sigma_y)*randn(1,K);

