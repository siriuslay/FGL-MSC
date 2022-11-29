function [W,W_n, W_obj, Z] = original_Update_W(mv, n, W, Z, K, Amv, Ai_Gi, alpha, beta, lambda, W_iter)
% 
for iter=1:W_iter
    W_obj = 0;
    parfor i=1:mv
        W(:,:,i)=round((K{i}+(lambda+beta)*eye(n)+eps)\(lambda*Z(:,:,i)+K{i}),4);
        T=W(:,:,i);
        T(find(T<0.00005))=0;
        T(logical(eye(size(T))))=zeros(1,n); %对角线变为0
        T = T./(sum(T,2)*ones(1,n)); %行归一化
        T=(T+T')/2;
        W(:,:,i)=T;
    end
    W_n = permute(W, [3,2,1]);
    Z_n = permute(Z, [3,2,1]);
    parfor ii=1:n
       Z_n(:,:,ii) = round((lambda*eye(mv)+alpha*Amv{ii}+eps)\(lambda*W_n(:,:,ii)+alpha*Ai_Gi{ii}),4); 
    end
    Z = permute(Z_n, [3,2,1]);
    parfor i = 1:mv
        T=Z(:,:,i);
        T(find(T<0.00005))=0;        
        T(logical(eye(size(T))))=zeros(1,n);
        T = T./(sum(T,2)*ones(1,n));
        T=(T+T')/2;
        Z(:,:,i)=T;
        W_obj = W_obj + (norm(Z(:,:,i)-W(:,:,i),'fro'));
    end
    if ((W_obj/mv)<1e-3)
        break
    end
    W=Z;
end
end
