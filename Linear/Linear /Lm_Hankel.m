% observe multiple nodes
close all
clear all
clc

N =10;
Create_J;
% find obervable nodes in all J1
for node = 1:N
    if all(NOT_OBSERVABLE_matrix(:,node) == 0)
        break;
    end
end


%%
[rank_J_N, rank_J_N_tol] = diyrank(J1^N);
for ai = 1:size(J,3)
    rank_estimation = [];
    estimation = [];
    ai
    J1 = J(:,:,ai);
    node_matrix = find(NOT_OBSERVABLE_matrix(ai,:)==0);
    for exp_t = 1:100
        if rank_J_N ~= N
            break;
        end
        observation = [];
        all_rank = [1];
        x = [];
        x(:,1)=1*unifrnd(-1,1,n,1); %this is x(0), with random i.i.d. entries unif(-5,5)
        
        A1_eig = eig(J1);
        %eps = 1e-15; % or eps=1e-10
        eps=1e-10;
        for t=1:3*n+200
            x(:,t+1) = J1 * x(:,t) + eps*normrnd(0,1,[n 1]) ;   %%version 1, with zero input
        end
        
        observation = x(node,1:end);
        max_matrix = [];
        for t=2:size(observation,2)/2-1
            % Hankel
            
            all_node_H = [];
            node_num = 5 ; % the number of observed nodes
            for node = 1:node_num
                H=[];
                observation = x(node_matrix(1,node),1:end);
                tmp_data = observation;
                for i1=1:t
                    for i2=1:t
                        H(i1,i2)=tmp_data(i1+i2-1);
                    end
                end
                all_node_H = [all_node_H,H];
            end
            
            % rank-func
            all_rank = [all_rank,diyrank(all_node_H)];
            
            % max-gap
            s = svd(all_node_H);
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
% for i = 1:100
%     for j = 1:100
%         ttt(i,j) = max(record_rank{i,j});
%     end
%     len_rank(i) = length(find(ttt(i,:)==N));
% end
% ttt = [];
% for i = 1:100
%     for j = 1:100
%         ttt(i,j) = max(record_max{i,j});
%     end
%     len_svd(i) = length(find(ttt(i,:)==N));
% end
Accuracy_rank = sum(len_rank)/(100*100)
Accuracy_max = sum(len_svd)/(100*100)