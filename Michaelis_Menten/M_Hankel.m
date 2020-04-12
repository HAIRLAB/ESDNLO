clear all
clc
close all

M_createdis;

node = 1; % select a node to measure

for exp_t = 1:S
    exp_t
    for ini = 1:100
        data = Z{exp_t,ini};
        observation = data(node,1:1:end);
        all_rank = [];
        max_matrix = [];
        
        for t=2:size(observation,2)/2-1 
            % Hankel
            H=[];
            sample_step = (floor(size(observation,2)/2/t));
            tmp_data = observation;
            
            for i1=1:t
                for i2=1:t
                    H(i1,i2)=tmp_data(i1+i2-1);
                end
            end
            
            all_rank = [all_rank rank(H)];
            
            s = svd(H);
            log_s = log(s);
            log_delta_s = log_s(1:end-1) - log_s(2:end);
            max_log_delta_s = find(log_delta_s==max(log_delta_s));
            max_matrix = [max_matrix max_log_delta_s];
        end
        rank_est = all_rank;
        
        Est_rank(exp_t,ini) = max(rank_est);
        Est_max(exp_t,ini) = max(max_matrix);
    end
end

Accuracy_rank = length(find(Est_rank==N))/(100*100)
Accuracy_svd = length(find(Est_max==N))/(100*100)