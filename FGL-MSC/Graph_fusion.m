function [Z_fu] = Graph_fusion(A, Z)
[n, mv] = size(A);
Z_fu = zeros(n);
for i = 1:mv
    Z_fu = Z_fu + diag(A(:,i))*Z(:,:,i);
end
end