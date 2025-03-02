function [Xmin,Xmax,dim,fobj] = fun_info(F)
switch F
    case 'F1'
       global model  
        N=5;%无人机的数量
        d=model.n;       % Number of Decision Variables = searching dimension of PSO = number of path nodes
        dim=3*d*N;%维度
        VarMax.r=2*norm(model.start-model.end)/d;           
        VarMin.r=0;
        % Inclination (elevation)
        AngleRange = pi/4; % Limit the angle range for better solutions
        VarMin.psi=-AngleRange;            
        VarMax.psi=AngleRange;          
        % Determine the angle of vector connecting the start and end points
        dirVector = model.end - model.start;
        phi0 = atan2(dirVector(2),dirVector(1));
        VarMin.phi=phi0 - AngleRange;           
        VarMax.phi=phi0 + AngleRange;
        Xmin=repmat([VarMin.r*ones(1,d),VarMin.psi*ones(1,d),VarMin.phi*ones(1,d)],[1,N]);%下限
        Xmax=repmat([VarMax.r*ones(1,d),VarMax.psi*ones(1,d),VarMax.phi*ones(1,d)],[1,N]);%上限
        fobj = @Cost;    % Cost Function; 
end
end
