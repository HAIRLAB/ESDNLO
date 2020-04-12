% Construct system matrix satisfying two assumptions
n = N;
J = zeros(N,N,100);
NOT_OBSERVABLE_matrix = zeros(100,N);

for j = 1:100 % 100 different system matrices
    step = 0.01;
    
    NOT_OBSERVABLE_FROM_NODE = zeros(1,N);
    diy_eigen_value = 0.98:-step:0.98-step*round(N/4); % make the convergence speed slow
    eigen_value = unifrnd (0.2,diy_eigen_value(end)-0.1 , N-1-length(diy_eigen_value),1);
    
    eigen_value = [1;diy_eigen_value';eigen_value];
    eigen_value = diag(eigen_value);
    
    right_eigen_vector = unifrnd(-1,1,N,N);
    right_eigen_vector(:,1) = ones(N,1);
    
    J1 = right_eigen_vector * eigen_value * inv(right_eigen_vector); % J1 satisfying two assumptions
    J(:,:,j) = J1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda = svd(J1);
    [xr, yr ] = eig(J1);
    
    v=inv(xr)';
    
    % test observable nodes of J1
    for u=1:n
        % C = e_u^T
        C=zeros(1,n);
        C(u)=1;
        % observability matrix O=[C;CJ;...CJ^{n-1}]
        O = C;
        for r=1:n-1
            O=[O;C*(J1^r)];
        end
        if rank(O) < n
            NOT_OBSERVABLE_FROM_NODE(1,u) = u;
        end
    end
    NOT_OBSERVABLE_matrix(j,:) = NOT_OBSERVABLE_FROM_NODE;  
end





