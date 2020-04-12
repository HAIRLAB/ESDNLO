function [output] = diynormalize(data)

[n,m]=size(data);
for i=1:n
    A(1,i)=norm(data(:,i));
end
A=repmat(A,n,1);
output=data./A;

end

