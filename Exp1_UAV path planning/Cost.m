function o= Cost(X)
global model
N=length(X)/(3*model.n);%���˻�������
o=0;
St=1;
for i=1:N %����ÿ�����˻�����Ӧ��ֵ
    Et=St+3*model.n-1;
    Xbest=X(St:Et);
    o=o+MyCost(Xbest,1);
    St=Et+1;
end
end