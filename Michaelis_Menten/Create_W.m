% Constructing W satisfying two assumptions
n = N;
W = zeros(N,N,100);
NOT_OBSERVABLE_matrix = zeros(100,N);
for j = 1:100
    step = 0.01;
    
    NOT_OBSERVABLE_FROM_NODE = zeros(1,N);
    diy_eigen_value = 0.98:-step:0.98-step*round(N/4);
    eigen_value = unifrnd (0.2,diy_eigen_value(end)-0.1 , N-1-length(diy_eigen_value),1);
    
    eigen_value = [1;diy_eigen_value';eigen_value];
    eigen_value = diag(eigen_value);
    
    right_eigen_vector = unifrnd(-1,1,N,N);
    right_eigen_vector(:,1) = ones(N,1);
    
    W1 = right_eigen_vector * eigen_value * inv(right_eigen_vector);
    W(:,:,j) = W1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda = svd(W1);
    [xr, yr ] = eig(W1);
    
    v=inv(xr)';
    for u=1:n
        %C = e_u^T
        C=zeros(1,n);
        C(u)=1;
        %observability matrix O=[C;CA;...CA^{n-1}]
        O = C;
        for r=1:n-1
            O=[O;C*(W1^r)];
        end
        if rank(O) < n
            NOT_OBSERVABLE_FROM_NODE(1,u) = u;
        end
    end
    NOT_OBSERVABLE_matrix(j,:) = NOT_OBSERVABLE_FROM_NODE;  
end





