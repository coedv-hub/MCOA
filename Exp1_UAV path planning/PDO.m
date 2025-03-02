
function [Best_PD,PDBest_P,PDConv]=PDO(N,T,LB,UB,Dim,F_obj)
PDBest_P=zeros(1,Dim);           % best positions
Best_PD=inf;                    % global best fitness
X=initializationPDO(N,Dim,UB,LB); %Initialize the positions of solution
Xnew=zeros(N,Dim);
PDConv=zeros(1,T);               % Convergance array
M=Dim;                      %set number of coteries 

t=1;                         % starting iteration
rho=0.005;                   % account for individual PD difference
% eps                         %food source quality 
epsPD=0.1;                  % food source alarm
OBest=zeros(1,size(X,1));     % old fitness values
CBest=zeros(1,size(X,1));     % new fitness values

for i=1:size(X,1) 
Flag_UB=X(i,:)>UB; % check if they exceed (up) the boundaries
Flag_LB=X(i,:)<LB; % check if they exceed (down) the boundaries
X(i,:)=(X(i,:).*(~(Flag_UB+Flag_LB)))+UB.*Flag_UB+LB.*Flag_LB;

    OBest(1,i)=F_obj(X(i,:));   %Calculate the fitness values of solutions
        if OBest(1,i)<Best_PD
            Best_PD=OBest(1,i);
            PDBest_P=X(i,:);
        end
end
  

while t<T+1  %Main loop %Update the Position of solutions
if mod(t,2)==0
             mu=-1;
else
             mu=1;
end
 
    DS=1.5*randn*(1-t/T)^(2*t/T)*mu;  % Digging strength
    PE=1.5*(1-t/T)^(2*t/T)*mu;  % Predator effect
    RL=levym(N,Dim,1.5);     % Levy random number vector
    TPD=repmat(PDBest_P,N,1); %Top PD
    for i=1:N 
        for j=1:M  
                cpd=rand*((TPD(i,j)-X(randi([1 N]),j)))/((TPD(i,j))+eps);
                P=rho+(X(i,j)-mean(X(i,:)))/(TPD(i,j)*(UB-LB)+eps);
                eCB=PDBest_P(1,j)*P;
                if (t<T/4)
                    Xnew(i,j)=PDBest_P(1,j)-eCB*epsPD-cpd*RL(i,j);    
                elseif (t<2*T/4 && t>=T/4)
                    Xnew(i,j)=PDBest_P(1,j)*X(randi([1 N]),j)*DS*RL(i,j);
                elseif (t<3*T/4 && t>=2*T/4)
                    Xnew(i,j)=PDBest_P(1,j)*PE*rand;
                else
                    Xnew(i,j)=PDBest_P(1,j)-eCB*eps-cpd*rand;
                end
        end
            
            Flag_UB=Xnew(i,:)>UB; % check if they exceed (up) the boundaries
            Flag_LB=Xnew(i,:)<LB; % check if they exceed (down) the boundaries
            Xnew(i,:)=(Xnew(i,:).*(~(Flag_UB+Flag_LB)))+UB.*Flag_UB+LB.*Flag_LB;
            CBest(1,i)=F_obj(Xnew(i,:));
            if CBest(1,i)<OBest(1,i)
                X(i,:)=Xnew(i,:);
                OBest(1,i)=CBest(1,i);
            end
            if OBest(1,i)<Best_PD
                Best_PD=OBest(1,i);
                PDBest_P=X(i,:);
            end
    end
  
    PDConv(t)=Best_PD;  %Update the convergence curve

    
     t=t+1;
end
end