function o= Cost(X)
global model
N=length(X)/(3*model.n);%无人机的数量
o=0;
St=1;
for i=1:N %计算每个无人机的适应度值
    Et=St+3*model.n-1;
    Xbest=X(St:Et);
    o=o+MyCost(Xbest,1);
    St=Et+1;
end
end