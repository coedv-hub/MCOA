close all
clear
clc
dbstop if all error
global model
model = CreateModel(); % 创建模型
F='F1';
[Xmin,Xmax,dim,fobj] = fun_info(F);%获取函数信息
pop=50;   %种群大小
maxgen=1;  %最大迭代次数
% [fMin1,bestX1,ConvergenceCurve1] = GRO(pop, maxgen,Xmin,Xmax,dim,fobj);
[fMin,bestX,ConvergenceCurve1] = MCOA(pop, maxgen,Xmin,Xmax,dim,fobj);

%% 计算无人机的相关信息
N=length(bestX)/(3*model.n);%无人机的数量
St=1;
for i=1:N %计算每个无人机的适应度值
    Et=St+3*model.n-1;
    Xbest=bestX(St:Et);
    BestPosition(i,:) = SphericalToCart(Xbest);%% 计算航迹坐标
    BestFit(i)=MyCost(Xbest,1);%% 计算每个无人机的适应度值
    UAVfit(i,:)=MyCost(Xbest,2);
    St=Et+1;
end

%% 保存结果
save BestPosition BestPosition %每个无人机的航迹坐标
save BestFit BestFit %每个无人机的总成本
save UAVfit UAVfit % 每个无人机的四个成本
save ConvergenceCurve ConvergenceCurve1 % 无人机集群的成本随迭代次数的变化

%% 画图
ColStr={'b-.','r--','c-.','m--','g-.'};%颜色
LegendStr={'UAV1','UAV2','UAV3','UAV4','UAV5'};

%图1 算法收敛曲线图
gca1=figure(1);
plot(ConvergenceCurve1,'r-','linewidth',3)
hold on
% plot(ConvergenceCurve2,'k','linewidth',3)
xlabel('迭代次数');
ylabel('全部无人机总成本');
legend({'GRO','GSEA'})

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


figure
bar(BestFit)
set(gca,'xtick',1:1:5);
set(gca,'XTickLabel',LegendStr)
ylabel('总成本')

figure
bar(UAVfit);
set(gca,'XTickLabel',LegendStr)
legend('路径成本','威胁成本','高度成本','转角成本')

figure
bar(UAVfit,"stacked");
set(gca,'XTickLabel',LegendStr)
legend('路径成本','威胁成本','高度成本','转角成本')

figure
bar(UAVfit');
set(gca,'xtick',1:1:4);
set(gca,'XTickLabel',{'路径成本','威胁成本','高度成本','转角成本'})
legend(LegendStr)

colormapStr=othercolor(61);
colormap(gca2,colormapStr);
colormap(gca3,colormapStr);

saveas(gca2,'Figure.fig');%将图二保存
openfig('Figure.fig');
view(2)
