
function drawPath(path,col)
%%%%
% xGrid=size(G,2);
% drawShanGe(G)
hold on
title('Grid Map Robot Path Planning Based on HEOA Algorithm','FontSize',12)
% set(gca,'XtickLabel',' ')
% set(gca,'YtickLabel',' ')
% xlabel('X坐标')
% ylabel('Y坐标')
L=size(path,1);
Sx=path(1,1)-0.5;
Sy=path(1,2)-0.5;
plot(Sx,Sy,'ro','MarkerSize',4,'LineWidth',4);   % 起点
for i=1:L-1
    figure(2)
    plot([path(i,2) path(i+1,2)]-0.5,[path(i,1) path(i+1,1)]-0.5,'Color',col,'LineWidth',2)
    hold on
%     %pause(0.5)
end
Ex=path(end,1)-0.5;
Ey=path(end,2)-0.5;
plot(Ey,Ex,'rs','MarkerSize',4,'LineWidth',4);   % 终点

