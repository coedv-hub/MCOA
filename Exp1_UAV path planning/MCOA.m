
function [best_fun,best_position,cuve_f,global_Cov]  =MCOA(N,T,lb,ub,dim,fobj)
%% Define Parameters
cuve_f=zeros(1,T); 
X=initialization(N,dim,ub,lb); %Initialize population
global_Cov = zeros(1,T);
Best_fitness = inf;
best_position = zeros(1,dim);
fitness_f = zeros(1,N);
for i=1:N
   fitness_f(i) =  fobj(X(i,:)); %Calculate the fitness value of the function
   if fitness_f(i)<Best_fitness
       Best_fitness = fitness_f(i);
       best_position = X(i,:);
   end
end
global_position = best_position; 
global_fitness = Best_fitness;
cuve_f(1)=Best_fitness;
t=1; 
while(t<=T)

    %% 加入折射对立-相互学习
    for i=1:N
        if rand<0.5
            Temp = X(i,:);
            alpha = max(Temp);
            beta = min(Temp);
            K = rand(1,dim);
            TempB = K.*(alpha + beta) - Temp; %求取反向解
            TempB = max(TempB, lb);
            TempB = min(TempB, ub);
            fitTemp = fobj(TempB);
            if(fitTemp<fitness_f(i))
                 X(i,:) = TempB;
                 fitness_f(i)=fitTemp;
            end     
        end
    end


    C = 2-(t/T); %Eq.(7)
    temp = rand*15+20; %Eq.(3)
    xf = (best_position+global_position)/2; % Eq.(5)
    Xfood = best_position;

    %% 重新定义洞穴
     Elite_pool(1,:)=best_position;
     Elite_pool(2,:)=global_position;
     Elite_pool(3,:)=(best_position+global_position)/2;
     Xmean=mean(X);
     Elite_pool(4,:)=Xmean;
     q = randi([10, N]);
     selected_index_q = randperm(N, q);
     Xq = X(selected_index_q, :);
     Xqmean = mean(Xq);
     Elite_pool(5,:)=Xqmean;
     

     %% 求两个中心
     % Randomly select 2 to 5 red-billed blue magpies
     p = randi([2, 5]);
     selected_index_p = randperm(N, p);
     Xp = X(selected_index_p, :);
     Xpmean = mean(Xp);
     Elite_pool(6,:)=Xpmean;

     %% 循环探索和开发
    for i = 1:N
        %% 探索
        if temp>30
            %% 在探索阶段加入其他的解，使其不易陷入局部最优
            if rand<0.5
                k1=randperm(6,1);
                Xnew(i,:) = X(i,:)+C*rand(1,dim).*(Elite_pool(k1,:)-X(i,:)); %Eq.(6)
            else
           

            %% 修改版的竞争洞穴
            for j=1:dim
            z = randperm(N,1);
            CF = (1 - t / T)^(2 * t / T);
            Xnew(i,j) = xf(j)+CF*(X(i,j)-X(z,j)).* randn(1, 1);
            end

               
            end
        else
            %% 开发
            %% foraging stage-觅食阶段
            P = 3*rand*(fitness_f(i)/fobj(Xfood)); %食物尺寸
            if P>2   % The food is too big  
                 Xfood = exp(-1/P).*Xfood;   %Eq.(12)
                for j = 1:dim
                    Xnew(i,j) = X(i,j)+(cos(2*pi*rand)-sin(2*pi*rand))*Xfood(j)*p_obj(temp); %Eq.(13)
                end
            else 
                Xnew(i,:) =((X(i,:)-Xfood))*p_obj(temp)+p_obj(temp).*rand(1,dim).*X(i,:); %Eq.(14)
            end
        end
    end
    %% boundary conditions
    for i=1:N
        for j =1:dim
            if length(ub)==1
                Xnew(i,j) = min(ub,Xnew(i,j));
                Xnew(i,j) = max(lb,Xnew(i,j));
            else
                Xnew(i,j) = min(ub(j),Xnew(i,j));
                Xnew(i,j) = max(lb(j),Xnew(i,j));
            end
        end
    end
   
    global_position = Xnew(1,:);
    global_fitness = fobj(global_position);
 
    for i =1:N
         %% Obtain the optimal solution for the updated population
        new_fitness = fobj(Xnew(i,:));
        if new_fitness<global_fitness
                 global_fitness = new_fitness;
                 global_position = Xnew(i,:);
        end
        %% Update the population to a new location
        if new_fitness<fitness_f(i)
             fitness_f(i) = new_fitness;
             X(i,:) = Xnew(i,:);
             if fitness_f(i)<Best_fitness
                 Best_fitness=fitness_f(i);
                 best_position = X(i,:);
             end
        end
    end
    global_Cov(t) = global_fitness;
    cuve_f(t) = Best_fitness;
    t=t+1;
%     if mod(t,50)==0
%       disp("COA"+"iter"+num2str(t)+": "+Best_fitness); 
%    end
end
 best_fun = Best_fitness;
end
function y = p_obj(x)   %Eq.(4)
    y = 0.2*(1/(sqrt(2*pi)*3))*exp(-(x-25).^2/(2*3.^2));
end