
clc
clear
close all
tic

outcomeMean=[];  % 保存均值
outcomebest=[];  % 保存最优值
outcomeworst=[];  % 保存最差值
%% 地图
G=[0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0; 
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
   0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0; 
   0 0 0 0 0 0 0 1 1 0 0 0 1 0 1 1 1 1 0 0; 
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0; 
   0 1 1 1 0 0 0 1 1 0 0 0 1 1 0 0 0 0 0 0; 
   0 1 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0;
   0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
   0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 0 0 0 0 0; 
   0 0 0 0 0 0 0 1 1 0 1 1 0 0 1 0 0 0 0 0; 
   0 0 0 0 0 0 0 1 1 0 1 1 0 0 1 0 0 0 0 0; 
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0; 
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0; 
   0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0; 
   0 0 0 1 0 0 1 1 0 0 0 1 0 0 0 0 0 0 0 0; 
   0 0 0 0 0 0 1 1 0 1 1 1 0 0 0 0 0 1 1 0; 
   0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0; 
   0 0 0 1 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0; 
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;];

num = size(G,1);
for i=1:num/2  
    for j=1:num
        m=G(i,j);
        n=G(num+1-i,j);
        G(i,j)=n;
        G(num+1-i,j)=m;
    end
end
%% 
S = [1 1];   
E = [num num];  
G0 = G;
G = G0(S(1):E(1),S(2):E(2)); 
[Xmax,dimensions] = size(G); X_min = 1;         
dimensions = dimensions - 2;            

%% 参数设置
max_gen = 1000;    % 最大迭代次数
num_polution = 50;         % 种群数量

fobj=@(x)fitness(x,G);


%% DBO
[Best_score1,bestX,DBO_curve]=DBO(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);
%% WOA
[Best_score2,~,WOA_curve]=WOA(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);
%% HLBO
% [Best_score3,~,HLBO_curve]=HLBO(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);
%% NRBO
[Best_score4,~,NRBO_curve]=NRBO(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);
%% LEA
[Best_score5,~,LEA_curve]=LEA(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);
%% HLOA
[Best_score6,~,HLOA_curve]=HLOA(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);
%% HEOA
[Best_score7,~,HEOA_curve]=HEOA(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);
%% GOOSE
[Best_score8,~,GOOSE_curve]=GOOSE(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);
%% CPO
[Best_score9,~,CPO_curve]=CPO(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);
%% GWO
[Best_score10,~,GWO_curve]=GWO(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);
%% BKA
[Best_score11,~,BKA_curve]=BKA(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);
%% COA
[Best_score12,~,COA_curve]=COA(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);
%% MCOA
[Best_score13,~,MCOA_curve]=MCOA(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);



disp(Best_score12)

