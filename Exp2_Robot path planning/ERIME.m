
function [Best_rime_rate,Best_rime,Convergence_curve]=ERIME(N,Max_iter,lb,ub,dim,fobj,G)

% initialize position
Best_rime=zeros(1,dim);
Best_rime_rate=inf;%change this to -inf for maximization problems

%% 初始化--改进1
Rimepop=initialization(N,dim,ub,lb);
% sine映射
% sine=3.8;
% Sine=rand(N,dim);
% for i=1:N
%     for j=2:dim
%         Sine(i,j)=(4/sine)*sin(pi*Sine(i,j-1));
%     end
% end
% result=Sine;
% Rimepop=repmat(lb,N,1)+result.*repmat((ub-lb),N,1);



Lb=lb.*ones(1,dim);% lower boundary 
Ub=ub.*ones(1,dim);% upper boundary
it=1;%Number of iterations
Convergence_curve=zeros(1,Max_iter);
Rime_rates=zeros(1,N);%Initialize the fitness value
newRime_rates=zeros(1,N);
W = 5;    %Soft-rime parameters, discussed in subsection 4.3.1 of the paper
%Initialization

%Calculate the fitness value of the initial position
for i=1:N
    Rime_rates(1,i)=fobj(Rimepop(i,:));
    %Calculate the fitness value for each search agent
    %Make greedy selections
    if Rime_rates(1,i)<Best_rime_rate
        Best_rime_rate=Rime_rates(1,i);
        Best_rime=Rimepop(i,:);
    end
end
for i = 1:N
    for j = 1:dim
       column = G(:,j+1);      % 地图的一列
       id = find(column == 0); % 该列自由栅格的位置
       x(i,j) =  id(randi(length(id))); % 随机选择一个自由栅格
       id = [];
    end
end
% Main loop
while it <= Max_iter
    RimeFactor = (rand-0.5)*2*cos((pi*it/(Max_iter/10)))*(1-round(it*W/Max_iter)/W);%Parameters of Eq.(3),(4),(5)
    E =(it/Max_iter)^0.5;%Eq.(6)
    newRimepop = Rimepop;%Recording new populations
    normalized_rime_rates=normr(Rime_rates);%Parameters of Eq.(7)
    %% 透镜成像
    for i=1:N
        Temp=Rimepop(i,:);
        k=(1+(it/Max_iter)^0.5)^10;
        TempB = (Ub+Lb)/2+(Ub+Lb)/(2*k)-Temp/k;
        TempB = max(TempB, lb);
        TempB = min(TempB, ub);
        fitTemp = fobj(TempB);
        if(fitTemp<Rime_rates(i))
             Rimepop(i,:) = TempB;
             Rime_rates(i)=fitTemp;
        end  
    end

    %% 计算质心
    for j=1:dim
        sum1=0;
        for i=1:N
            sum1=sum1+Rimepop(i,j);
        end
        X_centroid(j)=sum1/N;
    end

    for i=1:N
        for j=1:dim
            %Soft-rime search strategy
            r1=rand();
            if r1< E
               
                if rand > rand
                    newRimepop(i,j)=Best_rime(1,j)+RimeFactor*((Ub(j)-Lb(j))*rand+Lb(j));%Eq.(3)
                else
                    % XXX = max(Rimepop,[],1);
                    % MMM = min(Rimepop,[],1);
                    % newRimepop(i,j)=Best_rime(1,j)+RimeFactor*(((Ub(j)-Lb(j))*rand+Lb(j))+((XXX(j)-MMM(j))*rand+MMM(j)));%Eq.(3)
                    XXX = max(Rimepop,[],1);
                    MMM = min(Rimepop,[],1); %求每一维度的最小值
                    newRimepop(i,j)=Best_rime(1,j)+RimeFactor*((XXX(j)-MMM(j))*rand+MMM(j));%Eq.(3)
                end
            end
        end
    end

    %% 硬霜穿刺阶段
    for i=1:N
        for j=1:dim
            r2=rand;
            if r2 < normalized_rime_rates(i)
               r3=rand;
               RB=randn(N,dim);
               newRimepop(i,j)=Best_rime(1,j)+(r3*(Best_rime(1,j)-Rimepop(i,j))+(1-r3)*(X_centroid(j)-Rimepop(i,j)))*RB(i,j);%Eq.(3)
            end
        end
    end
    for i=1:N
        %Boundary absorption
        Flag4ub=newRimepop(i,:)>ub;
        Flag4lb=newRimepop(i,:)<lb;
        newRimepop(i,:)=(newRimepop(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        newRime_rates(1,i)=fobj(newRimepop(i,:));
        %Positive greedy selection mechanism
        if newRime_rates(1,i)<Rime_rates(1,i)
            Rime_rates(1,i) = newRime_rates(1,i);
            Rimepop(i,:) = newRimepop(i,:);
            if newRime_rates(1,i)< Best_rime_rate
               Best_rime_rate=Rime_rates(1,i);
               Best_rime=Rimepop(i,:);
            end
        end
    end
    Best_rime= LocalSearch(Best_rime,ub,G);
    Convergence_curve(it)=Best_rime_rate;
    it=it+1;
end





