clear all
clc
close all

K_createdis;

observable = [5,6]; % N=6
%observable = [1,5,9,10]; % N=10

for exp_t = 1:S
    data = THETA{exp_t};
    for index = 1:length(observable)
        node = observable(index);
        observation = data(node,1:1:end);
        all_rank = [];
        max_matrix = [];
        
        for t=2:size(observation,2)/2-1 % increase the size of Hankel
            % Hankel
            H=[];
            tmp_data = observation;
            
            for i1=1:t
                for i2=1:t
                    H(i1,i2)=tmp_data(i1+i2-1);
                end
            end
            
            % rank-func
            all_rank = [all_rank rank(H)];
            
            % max-gap
            s = svd(H);
            log_s = log(s);
            log_delta_s = log_s(1:end-1) - log_s(2:end) ;
            max_log_delta_s = find(log_delta_s==max(log_delta_s));
            max_matrix = [max_matrix max_log_delta_s];
        end
        rank_est = all_rank;
        
        Est_rank(exp_t,index) = max(rank_est);
        Est_max(exp_t,index) = max(max_matrix);
    end
end

Accuracy_rank = length(find(Est_rank==N))/(100*length(observable))
Accuracy_svd = length(find(Est_max==N))/(100*length(observable))

