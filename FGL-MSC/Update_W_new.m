function [W,W_n, W_obj, Z] = Update_W_new(mv, n, W, Z, X, A, G, alpha, beta, lambda, W_iter)
% 
for iter=1:W_iter
    W_obj = 0;
    for i=1:mv
        f=X{i};
        HH=f*f';
        W(:,:,i) = W(:,:,i).*( (2*HH + 2*lambda*Z(:,:,i))./(W(:,:,i)*HH + HH*W(:,:,i) + 2*lambda*W(:,:,i) + beta + eps) );
        % W矩阵的行归一化
        W(:,:,i) = W(:,:,i)./(sum(W(:,:,i),2)*ones(1,n));
        T=W(:,:,i);XX
        T(find(T<0))=0;
        T=(T+T')/2;
        T(logical(eye(size(T))))=zeros(1,n);
        W(:,:,i)=T;
    end
    W_n = permute(W, [3,2,1]);
    Z_n = permute(Z, [3,2,1]);
    for ii=1:n
       Amv = A(:,ii)*A(:,ii)';
       Ai_Gi = A(:,ii)*G(ii,:);
       Z_n(:,:,ii) = round((lambda*eye(mv)+alpha*Amv+eps)\(lambda*W_n(:,:,ii)+alpha*Ai_Gi),4); 
    end
    Z = permute(Z_n, [3,2,1]);
    for i = 1:mv
        T=Z(:,:,i);
        T(find(T<0))=0;
        T=(T+T')/2;
        T(logical(eye(size(T))))=zeros(1,n);
        Z(:,:,i)=T;
        W_obj = W_obj + (norm(Z(:,:,i)-W(:,:,i),'fro'));
    end
    if iter>5 && (W_obj<1e-3)
        break
    end
end

end