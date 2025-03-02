
function [best_score,best_pos,Convergence_curve]=GRO(N,Max_iter,lb,ub,dim,fobj )

lb = lb .* ones(1, dim);
ub = ub .* ones(1, dim);

%% GRO parameter initialization
sigma_initial = 2; 
sigma_final = 1 / Max_iter ;

% Initialize best position X* (global best)
best_pos=zeros(1, dim);
best_score=inf; %change this to -inf for maximization problems

%Initialize the gold prospectors? population Xi, i = 1, 2, . . . , N
Positions=initialization(N, dim, lb, ub);
Fit = inf(1,N);

%Initialize the gold prospectors? new positions Xnewi = Xi , i = 1, 2, . . . , N
X_NEW = Positions;
Fit_NEW = Fit;

Convergence_curve=zeros(1, Max_iter);
Convergence_curve(1) = min(Fit);
iter = 1;% Loop counter

%% Main loop
while iter <= Max_iter

    for i= 1:N   
        %Calculate fitness of current search agent at new position XNewi
        Fit_NEW(i) =  fobj(X_NEW(i,:));
        
        %Update position of current search agent Xi according to Equation (13)
        if Fit_NEW(i) < Fit(i)
            Fit(i) = Fit_NEW(i);
            Positions(i,:) = X_NEW(i,:);
        end
      
        %Update best search agent X*
        if Fit(i) < best_score
            % new gold mine is found
            best_score = Fit(i); 
            best_pos   = Positions(i,:);
        end

    end
    
   %Update l1, l2 by Equation (7)   
   l2 =  ((Max_iter - iter)/(Max_iter-1) )^2 * (sigma_initial - sigma_final) + sigma_final;
   l1 =  ((Max_iter - iter)/(Max_iter-1) )^1 * (sigma_initial - sigma_final) + sigma_final;
  
     
    %calculate the next position of current search agent XNewi with one of
    %... the migration, mining or collaboration methods
    for i = 1:size(Positions,1)

        coworkers = randperm(N-1,2);
        diggers = 1:N;
        diggers(i) = [];
        coworkers = diggers(coworkers);

        digger1 = coworkers(1);    %random prospector g1
        digger2 = coworkers(2);    %random prospector g2
        
        m = rand;
        %collaboration
        if m <  1/3
            for d  = 1:dim
                r1 = rand;                                         % r1 is a random number in [0,1]
                D3 = Positions(digger2,d) -  Positions(digger1,d); % Equation (11)
                X_NEW(i,d) = Positions(i,d) +  r1 * D3;            % Equation (12) 
            end
        %mining method
        elseif m < 2/3
            for d = 1:dim
                r1 = rand;                                          % r1 is a random number in [0,1]
                A2 = 2*l2*r1 - l2 ;                                 % Equation (10)               
                D2 = Positions(i,d) - Positions(digger1,d) ;        % Equation (8)
                X_NEW(i,d) = Positions(digger1,d) + A2*D2;          % Equation (9)                              
            end
        %migartion method    
        else
            for d = 1:dim
                r1 = rand; % r1 is a random number in [0,1]
                r2 = rand; % r2 is a random number in [0,1]
                C1 = 2 * r2;                                        % Equation (6)
                A1 = 1 + l1 * (r1 - 1/2);                           % Equation (5)
                D1 = C1 * best_pos(d) - Positions(i,d) ;            % Equation (3)
                X_NEW(i,d) = Positions(i,d) + A1 * D1;              % Equation (4)                      
            end
        end
            
        %Domain control
        X_NEW(i,:) = boundConstraint(X_NEW(i,:),Positions(i,:), lb , ub);
       
    end
    
    Convergence_curve(iter) = best_score;  
    display(['GRO:t=',num2str(iter),' fit=' num2str(best_score)]);
    iter = iter+1;
   
end


end


function newPos = boundConstraint(newPos, oldPos , lb, ub)

[NP, ~] = size(newPos);  % the Population size and the problem's dimension
%% check the lower bound
xl = repmat(lb, NP, 1);
pos = newPos < xl;
newPos(pos) =    oldPos(pos)  ;

%% check the upper bound
xu = repmat(ub, NP, 1);
pos = newPos > xu;
newPos(pos) =    oldPos(pos)  ;

end
function Positions=initialization(SearchAgents_no,dim,lb,ub)

% generate random position between lb and ub
delta = repmat(ub-lb,SearchAgents_no,1);
lb = repmat(lb, SearchAgents_no,1);
Positions = rand(SearchAgents_no,dim).* delta + lb;
end

