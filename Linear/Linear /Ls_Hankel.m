% observe a single node
close all
clear all
clc

N =10; % the number of nodes
Create_J;
% select an observable node
for node = 1:N
    if all(NOT_OBSERVABLE_matrix(:,node) == 0)
        break;
    end
end


%%
rank_J_N = rank(J1^N);

for ai = 1:size(J,3) % 100 different system matrices
    ai
    J1 = J(:,:,ai);
    rank_estimation = [];
    estimation = [];
    for exp_t = 1:100 % 100 different initial conditions
        if rank_J_N ~= N
            break;
        end
        observation = [];
        all_rank = [1];
        x = [];
        x(:,1)=1*unifrnd(-1,1,n,1); %this is x(0), with random i.i.d. entries unif(-1,1)
        
        eps = 1e-10; % noise intensity
        for t=1:3*n+200
            x(:,t+1) = J1 * x(:,t) + eps*normrnd(0,1,[n 1]) ; % version 1, with zero input
        end
        observation = x(node,1:end);
        max_matrix = [];
        for t=2:size(observation,2)/2-1 % increase the size of Hankel matrix
            
            % Hankel
            H=[];
            sample_step = (floor(size(observation,2)/2/t));
            tmp_data = observation;
            
            for i1=1:t
                for i2=1:t
                    H(i1,i2)=tmp_data(i1+i2-1);
                end
            end
            
            % rank-func
            all_rank = [all_rank,rank(H)];
            
            % max-gap
            s = svd(H);
            log_s = log(s);
            log_delta_s = log_s(1:end-1) - log_s(2:end) ;
            max_log_delta_s = find(log_delta_s==max(log_delta_s));
            max_matrix = [max_matrix max_log_delta_s];    
        end
        if  rank_J_N == N && max(all_rank)~= rank_J_N
            debug = 1;
        end
        record_rank{ai,exp_t} = all_rank;
        record_max{ai,exp_t} = max_matrix;
        rank_estimation = [rank_estimation max(all_rank)];
        estimation = [estimation max(max_matrix)];
        
    end
    % the number of correct estimation
    len_svd(ai) = length(find(estimation==N));
    len_rank(ai) = length(find(rank_estimation==N));
end

Accuracy_rank = sum(len_svd)/(100*100)
Accuracy_max = sum(len_rank)/(100*100)
