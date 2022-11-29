function [Z, Z_n] = Update_Z(Ai_Gi,Amv,W_n,alpha,lambda,n,mv)
    Z_n = W_n;
    parfor ii=1:n
       Z_n(:,:,ii) = (alpha*eye(mv)+lambda*Amv{ii}+eps)\(alpha*W_n(:,:,ii)+lambda*Ai_Gi{ii}); 
    end
    Z = permute(Z_n, [3,2,1]);
    parfor i = 1:mv
        T=Z(:,:,i);      
%         T(logical(eye(size(T))))=zeros(1,n);
        T = T./(sum(T,2)*ones(1,n));
        T(find(T<0.00005))=0;   
        T=(T+T')/2;
        Z(:,:,i)=T;
    end
    Z_n = permute(Z, [3,2,1]);
end