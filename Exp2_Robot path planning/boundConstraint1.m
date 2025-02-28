%This function is used for GRO bound checking
function newPos = boundConstraint1(newPos, oldPos , lb, ub)

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