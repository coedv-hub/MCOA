
clc
clear
close all
tic
%% 地图
G=[0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0; 
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
   1 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0; 
   1 1 0 1 0 0 1 1 0 0 0 1 0 0 0 0 0 0 0 0; 
   0 0 0 0 0 0 1 1 0 1 1 1 0 0 0 0 0 1 1 0; 
   0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0; 
   0 0 0 1 1 0 0 0 0 0 1 1 0 0 1 0 0 0 0 0; 
   0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0;];

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

%% GA
% [Best_score,Best_pos,GA_curve]=GA(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);
% %结果分析
% Best_pos = round(Best_pos);
% disp(['GA算法寻优得到的最短路径是：',num2str(Best_score)])
% route = [S(1) Best_pos E(1)];
% path_GA=generateContinuousRoute(route,G);
% path_GA=GenerateSmoothPath(path_GA,G);  
% path_GA=GenerateSmoothPath(path_GA,G);

%% SSA

% [Best_score,Best_pos,SSA_curve]=SSA(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);
% %结果分析
% Best_pos = round(Best_pos);
% disp(['SSA算法寻优得到的最短路径是：',num2str(Best_score)])
% route = [S(1) Best_pos E(1)];
% path_SSA=generateContinuousRoute(route,G);
% path_SSA=GenerateSmoothPath(path_SSA,G);  
% path_SSA=GenerateSmoothPath(path_SSA,G);

%% PSO

% [Best_score,Best_pos,PSO_curve]=PSO(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);
% %结果分析
% Best_pos = round(Best_pos);
% disp(['PSO算法寻优得到的最短路径是：',num2str(Best_score)])
% route = [S(1) Best_pos E(1)];
% path_PSO=generateContinuousRoute(route,G);
% path_PSO=GenerateSmoothPath(path_PSO,G);  
% path_PSO=GenerateSmoothPath(path_PSO,G);

%% DE

% [Best_score,Best_pos,DE_curve]=DE(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);
% %结果分析
% Best_pos = round(Best_pos);
% disp(['DE算法寻优得到的最短路径是：',num2str(Best_score)])
% route = [S(1) Best_pos E(1)];
% path_DE=generateContinuousRoute(route,G);
% path_DE=GenerateSmoothPath(path_DE,G);  
% path_DE=GenerateSmoothPath(path_DE,G);

%% GWO

% [Best_score,Best_pos,GWO_curve]=GWO(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);
% %结果分析
% Best_pos = round(Best_pos);
% disp(['GWO算法寻优得到的最短路径是：',num2str(Best_score)])
% route = [S(1) Best_pos E(1)];
% path_GWO=generateContinuousRoute(route,G);
% path_GWO=GenerateSmoothPath(path_GWO,G);  
% path_GWO=GenerateSmoothPath(path_GWO,G);

%% RIME

[Best_score,Best_pos,RIME_curve]=ERIME(num_polution,max_gen,X_min,Xmax,dimensions,fobj,G);
%结果分析
Best_pos = round(Best_pos);
disp(['RIME算法寻优得到的最短路径是：',num2str(Best_score)])
route = [S(1) Best_pos E(1)];
path_RIME=generateContinuousRoute(route,G);
path_RIME=GenerateSmoothPath(path_RIME,G);  
path_RIME=GenerateSmoothPath(path_RIME,G);


%% 画寻优曲线
% figure(1)
% plot(GA_curve,'k-o')
% hold on
% plot(SSA_curve,'y-^')
% hold on
% plot(PSO_curve,'b-*')
% hold on
% plot(DE_curve,'g-P')
% hold on
% plot(GWO_curve,'c-v')
% legend('GA','SSA','PSO','DE','GWO')
% title('20*20地图下各算法的收敛曲线')

%% 画路径
figure(2)
for i=1:num/2  
    for j=1:num
        m=G(i,j);
        n=G(num+1-i,j);
        G(i,j)=n;
        G(num+1-i,j)=m;
    end
end
 
n=num;
for i=1:num
    for j=1:num
        if G(i,j)==1 
            x1=j-1;y1=n-i; 
            x2=j;y2=n-i; 
            x3=j;y3=n-i+1; 
            x4=j-1;y4=n-i+1; 
             p=fill([x1,x2,x3,x4],[y1,y2,y3,y4],[0.2,0.2,0.2]); 
            % p.LineStyle="none";
            hold on 
        else 
            x1=j-1;y1=n-i; 
            x2=j;y2=n-i; 
            x3=j;y3=n-i+1; 
            x4=j-1;y4=n-i+1; 
            p=fill([x1,x2,x3,x4],[y1,y2,y3,y4],[1,1,1]); 
            % p.LineStyle="none";
            hold on 
        end 
    end 
end 
hold on

drawPath(path_RIME,'r')
% hold on
% drawPath(path_SSA,'y')
% hold on
% drawPath(path_PSO,'b')
% hold on
% drawPath(path_DE,'g')
% hold on
% drawPath(path_GWO,'c')



