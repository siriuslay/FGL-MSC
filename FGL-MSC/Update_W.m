function [W,W_n] = Update_W(mv, Z, W, K, alpha, beta)
    for iter = 1:5
        parfor i=1:mv
            W(:,:,i) = W(:,:,i).*( (2*K{i} + 2*alpha*Z(:,:,i))./(W(:,:,i)*K{i} + K{i}*W(:,:,i) + 2*alpha*W(:,:,i) + beta + eps) );
            T=W(:,:,i);
            T(find(T<0.00005))=0;
    %         T(logical(eye(size(T))))=zeros(1,n); %对角线变为0
    %         T = T./(sum(T,2)*ones(1,n)); %行归一化
            T=(T+T')/2;
            W(:,:,i)=T;
        end
    end
    W_n = permute(W, [3,2,1]);
end