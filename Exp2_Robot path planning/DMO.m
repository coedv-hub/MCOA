
function [BEF,BEP,BestCost]=DMO(nPop,MaxIt,VarMin,VarMax,nVar,F_obj,G)



%nVar=5;             % Number of Decision Variables

VarSize=[1 nVar];   % Decision Variables Matrix Size

%VarMin=-10;         % Decision Variables Lower Bound
%VarMax= 10;         % Decision Variables Upper Bound

%% ABC Settings

% MaxIt=1000;              % Maximum Number of Iterations

% nPop=100;               % Population Size (Family Size)

nBabysitter= 3;         % Number of babysitters

nAlphaGroup=nPop-nBabysitter;         % Number of Alpha group

nScout=nAlphaGroup;         % Number of Scouts

L=round(0.6*nVar*nBabysitter); % Babysitter Exchange Parameter 

peep=2;             % Alpha female鐥?vocalization 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Empty Mongoose Structure
empty_mongoose.Position=[];
empty_mongoose.Cost=[];

% Initialize Population Array
pop=repmat(empty_mongoose,nAlphaGroup,1);

% Initialize Best Solution Ever Found
BestSol.Cost=inf;
tau=inf;
Iter=1;
sm=inf(nAlphaGroup,1);

% Create Initial Population

for i=1:nAlphaGroup
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    for j=1:size(pop(i).Position,1)  
        Flag4ub=pop(i).Position>VarMax;
        Flag4lb=pop(i).Position<VarMin;
        pop(i).Position=(pop(i).Position.*(~(Flag4ub+Flag4lb)))+VarMax.*Flag4ub+VarMin.*Flag4lb; 
    end
       pop(i).Cost=F_obj(pop(i).Position);
    if pop(i).Cost<=BestSol.Cost
        BestSol=pop(i);
    end
end

% Abandonment Counter
C=zeros(nAlphaGroup,1);
CF=(1-Iter/MaxIt)^(2*Iter/MaxIt);

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);
for i = 1:nPop
    for j = 1:nVar
       column = G(:,j+1);      % 地图的一列
       id = find(column == 0); % 该列自由栅格的位置
       x(i,j) =  id(randi(length(id))); % 随机选择一个自由栅格
       id = [];
    end
end

%% DMOA Main Loop

for it=1:MaxIt
    
    % Alpha group
     F=zeros(nAlphaGroup,1);
     MeanCost = mean([pop.Cost]);
    for i=1:nAlphaGroup
        
        % Calculate Fitness Values and Selection of Alpha
        F(i) = exp(-pop(i).Cost/MeanCost); % Convert Cost to Fitness
    end
        P=F/sum(F);
      % Foraging led by Alpha female
    for m=1:nAlphaGroup
        
        % Select Alpha female
        i=RouletteWheelSelection(P);
        
        % Choose k randomly, not equal to Alpha
        K=[1:i-1 i+1:nAlphaGroup];
        k=K(randi([1 numel(K)]));
        
        % Define Vocalization Coeff.
        phi=(peep/2)*unifrnd(-1,+1,VarSize);
        % New Mongoose Position
        newpop.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
     for j=1:size(newpop.Position,1)  
        Flag4ub=newpop.Position>VarMax;
        Flag4lb=newpop.Position<VarMin;
        newpop.Position=(newpop.Position.*(~(Flag4ub+Flag4lb)))+VarMax.*Flag4ub+VarMin.*Flag4lb;
     end
        % Evaluation
        newpop.Cost=F_obj(newpop.Position);
        
        % Comparision
        if newpop.Cost<=pop(i).Cost
            pop(i)=newpop;
        else
            C(i)=C(i)+1;
        end
        
    end   
    
    % Scout group
    for i=1:nScout
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nAlphaGroup];
        k=K(randi([1 numel(K)]));
        
        % Define Vocalization Coeff.
        phi=(peep/2)*unifrnd(-1,+1,VarSize);
        
        % New Mongoose Position
        newpop.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
    for j=1:size(newpop.Position,1)  
        Flag4ub=newpop.Position>VarMax;
        Flag4lb=newpop.Position<VarMin;
        newpop.Position=(newpop.Position.*(~(Flag4ub+Flag4lb)))+VarMax.*Flag4ub+VarMin.*Flag4lb;
     end
        
        % Evaluation
        newpop.Cost=F_obj(newpop.Position);
        
        % Sleeping mould
        sm(i)=(newpop.Cost-pop(i).Cost)/max(newpop.Cost,pop(i).Cost);
        
        % Comparision
        if newpop.Cost<=pop(i).Cost
            pop(i)=newpop;
        else
            C(i)=C(i)+1;
        end
        
    end    
    % Babysitters
    for i=1:nBabysitter
        if C(i)>=L
            pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
       for j=1:size(pop(i).Position,1)  
        Flag4ub=pop(i).Position>VarMax;
        Flag4lb=pop(i).Position<VarMin;
        pop(i).Position=(pop(i).Position.*(~(Flag4ub+Flag4lb)))+VarMax.*Flag4ub+VarMin.*Flag4lb;
       end
            pop(i).Cost=F_obj(pop(i).Position);
            C(i)=0;
        end
    end    
     % Update Best Solution Ever Found
    for i=1:nAlphaGroup
        if pop(i).Cost<=BestSol.Cost
            BestSol=pop(i);
        end
    end    
        
   % Next Mongoose Position
   newtau=mean(sm);
   for i=1:nScout
        M=(pop(i).Position.*sm(i))/pop(i).Position;
        if newtau>tau
           newpop.Position=pop(i).Position-CF*phi*rand.*(pop(i).Position-M);
        else
           newpop.Position=pop(i).Position+CF*phi*rand.*(pop(i).Position-M);
        end
        tau=newtau;
      for j=1:size(newpop.Position,1)  
        Flag4ub=newpop.Position>VarMax;
        Flag4lb=newpop.Position<VarMin;
        newpop.Position=(newpop.Position.*(~(Flag4ub+Flag4lb)))+VarMax.*Flag4ub+VarMin.*Flag4lb;
      end
   end
       
   % Update Best Solution Ever Found
    for i=1:nAlphaGroup
        if pop(i).Cost<=BestSol.Cost
            BestSol=pop(i);
        end
    end
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    BEF=BestSol.Cost;
    BEP=BestSol.Position;
    BEP= LocalSearch(BEP,VarMax,G);
    % Display Iteration Information
%     disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    

end
end

function i=RouletteWheelSelection(P)

    r=rand;
    
    C=cumsum(P);
    
    i=find(r<=C,1,'first');

end


