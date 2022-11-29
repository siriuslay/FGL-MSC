function [A] = Update_A(A, G, Z_n)
[n, mv] = size(A);
parfor i=1:n
   H = round(G(:,i)*ones(1,mv)-Z_n(:,:,i)',4);
   HH = H'*H;
   HH = HH + eye(mv)*eps*trace(HH); %保证HH可逆
   T = HH\ones(mv,1);
   T(find(T<0))=0;
   A(i,:) = T/sum(T);
%    A(i,:)=round(((HH*ones(mv,1))/(ones(1,mv)*HH*ones(mv,1)))',4);
end
end