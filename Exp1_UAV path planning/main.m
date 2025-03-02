close all
clear
clc
dbstop if all error
global model
model = CreateModel(); % ����ģ��
F='F1';
[Xmin,Xmax,dim,fobj] = fun_info(F);%��ȡ������Ϣ
pop=50;   %��Ⱥ��С
maxgen=1;  %����������
% [fMin1,bestX1,ConvergenceCurve1] = GRO(pop, maxgen,Xmin,Xmax,dim,fobj);
[fMin,bestX,ConvergenceCurve1] = MCOA(pop, maxgen,Xmin,Xmax,dim,fobj);

%% �������˻��������Ϣ
N=length(bestX)/(3*model.n);%���˻�������
St=1;
for i=1:N %����ÿ�����˻�����Ӧ��ֵ
    Et=St+3*model.n-1;
    Xbest=bestX(St:Et);
    BestPosition(i,:) = SphericalToCart(Xbest);%% ���㺽������
    BestFit(i)=MyCost(Xbest,1);%% ����ÿ�����˻�����Ӧ��ֵ
    UAVfit(i,:)=MyCost(Xbest,2);
    St=Et+1;
end

%% ������
save BestPosition BestPosition %ÿ�����˻��ĺ�������
save BestFit BestFit %ÿ�����˻����ܳɱ�
save UAVfit UAVfit % ÿ�����˻����ĸ��ɱ�
save ConvergenceCurve ConvergenceCurve1 % ���˻���Ⱥ�ĳɱ�����������ı仯

%% ��ͼ
ColStr={'b-.','r--','c-.','m--','g-.'};%��ɫ
LegendStr={'UAV1','UAV2','UAV3','UAV4','UAV5'};

%ͼ1 �㷨��������ͼ
gca1=figure(1);
plot(ConvergenceCurve1,'r-','linewidth',3)
hold on
% plot(ConvergenceCurve2,'k','linewidth',3)
xlabel('��������');
ylabel('ȫ�����˻��ܳɱ�');
legend({'GRO','GSEA'})

%ͼ2��ͼ3 ���˻��켣ͼ
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
ylabel('�ܳɱ�')

figure
bar(UAVfit);
set(gca,'XTickLabel',LegendStr)
legend('·���ɱ�','��в�ɱ�','�߶ȳɱ�','ת�ǳɱ�')

figure
bar(UAVfit,"stacked");
set(gca,'XTickLabel',LegendStr)
legend('·���ɱ�','��в�ɱ�','�߶ȳɱ�','ת�ǳɱ�')

figure
bar(UAVfit');
set(gca,'xtick',1:1:4);
set(gca,'XTickLabel',{'·���ɱ�','��в�ɱ�','�߶ȳɱ�','ת�ǳɱ�'})
legend(LegendStr)

colormapStr=othercolor(61);
colormap(gca2,colormapStr);
colormap(gca3,colormapStr);

saveas(gca2,'Figure.fig');%��ͼ������
openfig('Figure.fig');
view(2)
