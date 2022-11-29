clear;clc;
% f=load('reuters_1200.mat');%=load('C:\Users\User\Desktop\research\multiviewlearning\shiguoxin\bbc_seg14of4.mat');
% data=f.data;
% label=f.labels;

f=load('MSRC-v1.mat');
data=f.X';
label=f.Y + 1;

% f=load('Caltech101-7.mat');
% data=f.X';
% label=f.Y;

% f=load('ORL_mtv.mat');  %10 1000 0.0001 10 66.525
% label=f.gt;
% for i = 1:size(f.X,2)
%     data{i}=f.X{i}';
% end

% f=load('miniHW.mat');  % 
% data=f.X';
% label=f.Y;

% f=load('handwritten.mat');  % 
% data=f.X';
% label=f.Y+1;

% f=load('bbcsport.mat');  % 1 12 0.001(10)-85.974
% data=f.X';
% label=f.Y;

% f=load('new_WebKB.mat');  % 1000	1000	0.01	10 -89.72
% data=f.data;
% label=f.label;

n=size(data{1},1);
beta = 1;
lambda = 1;
k = 5;

parfor i=1:size(data,1)
dist = max(max(data{i})) - min(min(data{i}));
m01 = (data{i} - min(min(data{i})))/dist;
data{i} = 2 * m01 - 1;
K{i}=data{i}*data{i}';

% Z{i} = constructW_PKN(data{i}', k, 1); %knn
W{i} = (K{i}+ (beta+lambda)* eye(n)) \ (K{i} + lambda*eye(n));
end



% 
% 
% c=length(unique(true_label));
% G = zeros(n);
% W=repmat(ones(n)-eye(n),[1,1,mv]);
% A = ones(n,mv)/mv;


