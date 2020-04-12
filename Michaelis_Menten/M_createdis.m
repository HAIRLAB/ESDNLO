clear all
clc

%% topograph
n = 10;
N = 10;
Create_W;
W1 = W;
%% model parameters
N = n;
S = 100;
epsilon = 1/N;

%% system evolution
for exp_t = 1:S
    W = W1(:,:,exp_t);
    for ini = 1:100
        z = [];
        init = unifrnd(-0.00000000000001,0.00000000000001,n,1);
        z(:,1) = init;
        for k = 1:3*n+400
            for i = 1:n
                tmp = 0;
                for j = 1:n
                    tmp = tmp + W(i,j) * z(j,k)/(1+z(j,k));
                end
                z(i,k+1) = (1-epsilon) * z(i,k) + epsilon * tmp;
            end
        end
        Z{exp_t,ini} = z;
    end
end

