
clear;clc;

f=load('data\Caltech20.mat');   % 100leaves: 0.01	0.01	0.1	10    | ORL: 0.01 1 1 10
data=f.X;
label=f.Y;

para1=[0.01 0.1 1 ];
para3=[0.1 1 10];
para4=[0.01 0.1 1 10];

parfor i=1:size(data,1)
dist = max(max(data{i})) - min(min(data{i}));
data{i} = (data{i} - min(min(data{i})))/dist;
% data{i} = 2 * m01 - 1;
K{i}=data{i}*data{i}';
end

% [result,W,G,A,F, F_diff, obj] = Noisy_fu(data,K,label, 100000000000, 1, 1, 10);
% dlmwrite('results\ablation_ORL.txt',[100000000000 1 1 10 result(7,:) result(4,:) result(1,:) result(5,:) F_diff obj],'-append','delimiter','\t','newline','pc');
% [result,W,G,A,F, F_diff, obj] = FGL(data,K,label, 0.01, 1, 1, 0);
% dlmwrite('results\ablation_ORL.txt',[0.01 1 1 0 result(7,:) result(4,:) result(1,:) result(5,:) F_diff obj],'-append','delimiter','\t','newline','pc');
[result,W,G,A,F, F_diff, obj] = FMSC(data,K,label,0.1, 0.01, 0.1, 10);
dlmwrite('results\FMSC_Caltech_conv.txt',[0.1 0.01 0.1 10 result(7,:) result(4,:) result(1,:) result(5,:) F_diff obj],'-append','delimiter','\t','newline','pc');

% for i=1:length(para1)
%     for j=i
%         for k=1:length(para3)
%             for t=1:length(para4)
%                 [result,W,G,A,F] = FMSC(data,K,label,para1(i),para1(j),para3(k),para4(t));
%                 dlmwrite('results\FMSC_ORL_conv.txt',[para1(i) para1(j) para3(k) para4(t) result(7,:)],'-append','delimiter','\t','newline','pc');
%             end
%         end
%     end
% end

clear;clc;