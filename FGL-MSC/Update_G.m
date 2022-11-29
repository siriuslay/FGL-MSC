function [G] = Update_G(Z_fu, F, lambda, gamma, n,k)
dist_F = L2_distance_1(F',F');
G=zeros(n);

D = (gamma/2)*dist_F-2*lambda*Z_fu;
[dumb, idx] = sort(D, 2); % sort each row

for i = 1:n
    id = idx(i,2:k+2);
    di = D(i, id);
    G(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end
end