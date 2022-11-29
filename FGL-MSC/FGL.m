function [result,W,G,A,F, F_diff, obj] = FGL(X,K,true_label,alpha,beta,lambda,gamma)
%FMSC 细粒度多图融合谱聚类
% sum||Xv-Xv*Wv||_F + alpha*||Wv-Zv||_F + beta*||W||_1 + lambda*||G-sum(Av*Zv)||_F +η*||G||_F + gamma*Tr(F'LF)

% W:①||Xv-Xv*Wv||_F + alpha*||Wv-Zv||_F + beta*||W||_1
% Z:②alpha*||Wv-Zv||_F + lambda*||G-sum(Av*Zv)||_F
% G:③lambda*||G-sum(Av*Zv)||_F + gamma*Tr(F'LF)
% A:④||G-sum(Av*Zv)||_F
% F:⑤Tr(F'LF)
%   X: a cell of feature matrix,n*d
%   true_label: the true labels of samples
%   alpha,beta,gamma: hyperparameters
%   A:n*mv F:n*c

mv=size(X,1);
n=size(X{1},1);
c=length(unique(true_label));
G = eye(n);
Z = repmat(eye(n),[1,1,mv]);
A = ones(n,mv)/mv;
% ===============
W = repmat(ones(n)-eye(n),[1,1,mv]);     %= ==============

maxIter = 10;
tol = 1e-3;
k=10;

for i = 1:mv
    Z(:,:,i) = constructW_PKN(X{i}', k, 1); %knn
    G = G + diag(A(:,i))*Z(:,:,i);
end
G = (G+G')/2;
D = diag(sum(G));
L = D-G; 
[F, temp, ev]=eig1(L, c, 0);
F_diff = zeros(1, maxIter);
obj = zeros(1,maxIter);

for iter = 1:maxIter
    fprintf('Iteration: %d\n', iter);
    F_old = F;
    parfor i = 1:n
       Amv{i} = A(i,:)'*A(i,:);
       Ai_Gi{i} = A(i,:)'*G(i,:);
    end
    % Update W
    [W,W_n] = Update_W(mv, Z, W, K, alpha, beta);    
    % Update Z
    [Z, Z_n] = Update_Z(Ai_Gi,Amv,W_n,alpha,lambda,n,mv);
    % multi-graph fusion
    [Z_fu] = Graph_fusion(A, Z);
    % Update G
    [G] = Update_G(Z_fu, F, lambda, gamma, n,k);
    G = (G+G')/2;
    D = diag(sum(G));
    L = D-G;
    % Update A
    [A] = Update_A(A, G, Z_n);
    % Update F
%     [F, temp, ev]=eig1(L, c, 0);
    
    % output
    loss1 = 0;loss2 = 0;
    parfor i = 1:mv
        loss1 = loss1 + norm((X{i}-W(:,:,i)*X{i}), 'fro') + alpha* norm((W(:,:,i)-Z(:,:,i)), 'fro') + beta*norm(W(:,:,i), 1);
    end
    loss2 = lambda * norm(G-Z_fu,"fro") + gamma * norm(G, "fro");  
    obj(1,iter) = loss1 + loss2;
    F_diff(1,iter) = sum(sum(abs(F-F_old)));
    if F_diff(1,iter) < tol
        break;
    end
end

[F, temp, ev]=eig1(L, c, 0);

res=zeros(10,8);
for ij=1:10
    actual_ids= kmeans(F, c, 'emptyaction', 'singleton', 'replicates', 1, 'display', 'off');
    [res(ij,:)] = Clustering8Measure( actual_ids,true_label);
end
% result(1,1)='Fscore';
% result(2,1)='Precision';
% result(3,1)='Recall';
% result(4,1)='nmi';
% result(5,1)='AR';
% result(6,1)='Entropy';
% result(7,1)='ACC';
% result(8,1)='Purity';
result(1,1)=mean(res(:,1));result(1,2)=std(res(:,1));
result(2,1)=mean(res(:,2));result(2,2)=std(res(:,2));
result(3,1)=mean(res(:,3));result(3,2)=std(res(:,3));
result(4,1)=mean(res(:,4));result(4,2)=std(res(:,4));
result(5,1)=mean(res(:,5));result(5,2)=std(res(:,5));
result(6,1)=mean(res(:,6));result(6,2)=std(res(:,6));
result(7,1)=mean(res(:,7));result(7,2)=std(res(:,7));
result(8,1)=mean(res(:,8));result(8,2)=std(res(:,8));