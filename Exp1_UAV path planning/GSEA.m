function [best_Score,best_pos,convergence_curve]=GSEA(N,Tmax,lb,ub,dim,fobj)


%% 初始化最优适应度值，最优解以及收敛曲线
best_pos=[];
best_Score=inf;
convergence_curve=zeros(1,Tmax);
%% 种群初始化--使用随机初始化
% 第一种映射方式--随机初始化
X=initialization(N,dim,ub,lb);
% 计算适应度，更新最优解
for i=1:N
    pos_fit(i)=fobj(X(i,:));%#ok
    if pos_fit(i)<best_Score 
            best_Score=pos_fit(i); 
            best_pos=X(i,:);
    end
end

%% 定义一个Xnew，用于进行位置更新
Xnew=X(:,:);
%% 定义上下界
lb=ones(1,dim).*lb; 
ub=ones(1,dim).*ub; 

%% 初始化控制参数
Elite_pool=[]; 

%% 主循环
for t=1:Tmax 
    %% 探索阶段
    % 确定精英池
    [~, idx1]=sort(X);
    best_pos=X(idx1(1),:);
    second_pos=X(idx1(2),:);
    third_pos=X(idx1(3),:);
    fourth_pos=X(idx1(4),:);
    Xmean=mean(X);
    P1=N/2+10;  
    P2=N/2;  % 用于保存中位数位置
    Worst_pos=X(idx1(P1),:);
    median_pos=X(idx1(P2),:); 
    Elite_pool(1,:)=second_pos;
    Elite_pool(2,:)=third_pos;
    Elite_pool(3,:)=fourth_pos;
    Elite_pool(4,:)=Xmean;
    Elite_pool(5,:)=median_pos;
    Elite_pool(6,:)=Worst_pos;
    % 探索阶段参数设置
    
    %探索公式
    for i=1:N
        RB=randn(N,dim);          %Brownian random number vector
        r1=rand;  % 随机因子
        r2=rand;  % 随机因子
        k1=randperm(6,1);  
        k2=randperm(6,1);  
        T1=randperm(N,1);  
 
        Xnew(i,:)=X(i,:)+(r1*(Elite_pool(k1,:)-X(i,:))+(1-r1-r2)*(X(T1,:)-X(i,:))+r2*(best_pos-X(i,:))).* randn(1, dim);    
    end
    %边界控制
    Xnew = boundaryCheck(Xnew, lb, ub);
    %种群更新
    for i=1:N
         New_fit= fobj(Xnew(i,:));
         if New_fit<pos_fit(i)
            pos_fit(i)=New_fit;
            X(i,:)=Xnew(i,:);
         end
         if New_fit<best_Score
            best_Score=New_fit; 
            best_pos=Xnew(i,:);
         end
    end

    %% 开发阶段
    N1=N*0.8;  
    N2=N*0.2;  
    for i=1:N
         index(i)=i;
    end
    index1=randperm(N,N1);
    index2=setdiff(index,index1);
    % 开发阶段参数
    TF = 1+(t/Tmax); % levy飞行动态参数（保留做备用）
   
    % 螺旋飞行函数
     k=1; 
     L = 2*rand-1;
     z = exp(k*cos(pi*(1-(t/Tmax))));
     %% 开发阶段位置更新
    for i=1:N1
        if rand<rand  
            A1=randperm(N,1);  
            T2=randperm(N1,1);  
            k3=randperm(6,1); 
            GZ=N1/3;
            index3=randperm(N1,GZ);
            B1=randperm(GZ,1);  
            Xnew(index1(i),:)=X(i,:)+(((best_pos-X(i,:))+(Elite_pool(k3,:)-X(i,:))+(X(index3(B1),:)-X(i,:))+(X(A1,:)-X(i,:))+(X(index1(T2),:)-X(i,:)))/5)*exp(z*L)*cos(2*pi*L);
        else  
            T3=randperm(N1,1);  
            k4=randperm(6,1); 
            GZ=N1/3;
            index3=randperm(N1,GZ);
            B2=randperm(GZ,1); 
            Xnew(index1(i),:)=X(i,:)+(((best_pos-X(i,:))+(Elite_pool(k4,:)-X(i,:))+(X(index3(B2),:)-X(i,:))+(X(index1(T3),:)-X(i,:)))/4).*Levy(dim)*TF;
        end
    end
  
    if N1<N
        N1=N1+1;
        N2=N2-1;
    end
   
    if N2>=1
        for i=1:N2
            if rand <rand  
                A2=randperm(N,1);  
                k5=randperm(6,1); 
                k6=randperm(N2,1);
                T4=randperm(N,1); 
                Xnew(index2(i),:)=X(i,:)+(((best_pos-X(i,:))+(Elite_pool(k5,:)-X(i,:))+(X(index2(k6),:)-X(i,:))+(X(A2,:)-X(i,:))+(X(T4,:)-X(i,:)))/5)*exp(z*L)*cos(2*pi*L);
            else  
                k7=randperm(6,1); 
                k8=randperm(N2,1); 
                T5=randperm(N,1);  
                Xnew(index2(i),:)=X(i,:)+(((best_pos-X(i,:))+(Elite_pool(k7,:)-X(i,:))+(X(index2(k8),:)-X(i,:))+(X(T5,:)-X(i,:)))/4).*Levy(dim)*TF;
            end
        end
    end
    
    %% 边界控制
    Xnew = boundaryCheck(Xnew, lb, ub);
    %种群更新
    for i=1:N
         New_fit= fobj(Xnew(i,:));
         if New_fit<pos_fit(i)
            pos_fit(i)=New_fit;
            X(i,:)=Xnew(i,:);
         end
         if New_fit<best_Score
            best_Score=New_fit; 
            best_pos=Xnew(i,:);
         end
    end

    %% 贪婪选择--防止种群退化
    for i=1:N
        newfit(1,i)=fobj(Xnew(i,:));
        if newfit(1,i)<pos_fit(1,i)
            pos_fit(1,i) = newfit(1,i);
            X(i,:) = Xnew(i,:);
            if newfit(1,i)< best_Score
               best_Score=newfit(1,i);
               best_pos=Xnew(i,:);
            end
        end
    end
    convergence_curve(t)=best_Score;
end

end







%%  随机初始化函数
function [ X ]=initialization(N,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    X=rand(N,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        X(:,i)=rand(N,1).*(ub_i-lb_i)+lb_i;
    end
end
end

%% 边界控制函数
function [ X ] = boundaryCheck(X, lb, ub)

    for i=1:size(X,1)
            FU=X(i,:)>ub;
            FL=X(i,:)<lb;
            X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
    end
end

%% levy 飞行函数
function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end




