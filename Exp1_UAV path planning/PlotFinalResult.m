close all
clear
clc
warning ('off')
global model
model = CreateModel(); % Create search map and parameters
load('BestPosition.mat');
load('BestFit.mat');
load('ConvergenceCurve.mat');
load('UAVfit.mat');

%% 画图
ColStr={'b-.','r--','c-.','m--','g-.'};%算法代码的颜色
LegendStr={'UAV1','UAV2','UAV3','UAV4','UAV5'};

%图1 算法收敛曲线图
gca1=figure(1);
plot(ConvergenceCurve1,'r-','linewidth',3)
xlabel('迭代次数');
ylabel('全部无人机总成本');
legend('GRO')

%图2和图3 无人机轨迹图
gca2=PlotModel(model);
gca3=figure(3);


[h11,h12]=PlotSolution(BestPosition(1,:),model,ColStr{1},gca2,gca3);
[h21,h22]=PlotSolution(BestPosition(2,:),model,ColStr{2},gca2,gca3);
[h31,h32]=PlotSolution(BestPosition(3,:),model,ColStr{3},gca2,gca3);
[h41,h42]=PlotSolution(BestPosition(4,:),model,ColStr{4},gca2,gca3);
[h51,h52]=PlotSolution(BestPosition(5,:),model,ColStr{5},gca2,gca3);
legend([h11,h21,h31,h41,h51],LegendStr,'location','NorthWest');
legend([h12,h22,h32,h42,h52],LegendStr,'location','NorthWest');




