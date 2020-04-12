% create trajectories of coupled Kuramoto oscilators
clear all
clc

%% topograph (N=6)
% the 5th and 6th nodes are observable 
n = 6;
A = [0,1,0,1,1,1;
    1,0,1,0,0,0;
    0,1,0,0,0,0;
    1,0,0,0,0,0;
    1,0,0,0,0,1;
    1,0,0,0,1,0]; % adjaceny matrix
%% topograph (N=10)
% the 1st, 9th and 10th nodes are observable 
% n = 10;
% A = [0,1,0,0,0,0,1,0,0,0;
%      1,0,1,0,1,1,0,1,0,0;
%      0,1,0,1,0,0,0,0,0,0;
%      0,0,1,0,0,0,0,1,1,1;
%      0,1,0,0,0,1,0,0,0,0;
%      0,1,0,0,1,0,0,0,1,0;
%      1,0,0,0,0,0,0,1,0,1;
%      0,1,0,1,0,0,1,0,0,0;
%      0,0,0,1,0,1,0,0,0,0;
%      0,0,0,1,0,0,1,0,0,0];
%% model parameters
N = n;
K = N;
w = unifrnd(-0.0000000001,0.0000000001,N-1,1);
w = [w; -sum(w)];
epsilon = 1/N;
S = 100;

%% system evolution
for exp_t = 1:S
    theta = [];
    init = unifrnd(-0.000001,0.000001,N,1);
    theta(:,1) = init;
    for k = 1:3*n+200
        for i = 1:n
            tmp = 0;
            for j = 1:n
                tmp = tmp + A(i,j) * sin(theta(j,k) - theta(i,k));
            end
            theta(i,k+1) = theta(i,k) + epsilon * (w(i) + K/N * tmp);
        end
    end
    THETA{exp_t} = theta;
end

