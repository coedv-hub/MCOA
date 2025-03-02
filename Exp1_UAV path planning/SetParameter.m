%% 参数设置
% (1) 设置每个无人机的待优化点的个数   
n=8;
% (2) 设置起点和终点
start_location = [50;150;50];
end_location = [900;720;150];

% (3)设置障碍物 
R1=30;  % Radius 80
x1 = 400; y1 = 500; z1 = 150; % center

R2=50;  % Radius 70
x2 = 700; y2 = 150; z2 = 150; % center

R3=40;  % Radius 80
x3 = 550; y3 = 450; z3 = 150; % center

R4=50;  % Radius 70
x4 = 350; y4 = 100; z4 = 150; % center

R5=30;  % Radius 70
x5 = 400; y5 = 650; z5 = 150; % center

R6=30;  % Radius 80
x6 = 800; y6 = 800; z6 = 150; % center

R7=70;  % Radius 70
x7 = 750; y7 = 350; z7 = 150; % center

R8=60; 
x8 = 150; y8 = 350; z8 = 150; % center

R9=50; 
x9 = 920; y9 = 600; z9 = 150; % center

R10=50; 
x10 = 920; y10 = 200; z10 = 150; % center
Threats=[x1 y1 z1 R1;x2 y2 z2 R2; x3 y3 z3 R3; x4 y4 z4 R4; x5 y5 z5 R5;x6 y6 z6 R6;x7 y7 z7 R7;x8 y8 z8 R8;x9 y9 z9 R9;x10 y10 z10 R10];