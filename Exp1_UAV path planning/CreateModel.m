%{
This model is generated by:
- Loading terrain map
- Creating threats as cylinders
- Creating start and finish points
- Setting ranges and limits
%}

function model=CreateModel()
    SetParameter;%参数设置
%     H = imread('ChrismasTerrain.tif'); % Get elevation data
%     H (H < 0) = 0;
    [~,~,z]=peaks(1000);%产生山峰地图
    H=30*abs(z);
    
    MAPSIZE_X = size(H,2); % x index: columns of H
    MAPSIZE_Y = size(H,1); % y index: rows of H
    [X,Y] = meshgrid(1:MAPSIZE_X,1:MAPSIZE_Y); % Create all (x,y) points to plot


    
    % Map limits
    xmin= 1;
    xmax= MAPSIZE_X;
    
    ymin= 1;
    ymax= MAPSIZE_Y;
    
    zmin = 100;
    zmax = 200;  
  
    % Number of path nodes (not including the start position (start node))

    
    % Incorporate map and searching parameters to a model
    model.start=start_location;
    model.end=end_location;
    model.n=n;
    model.xmin=xmin;
    model.xmax=xmax;
    model.zmin=zmin;
    model.ymin=ymin;
    model.ymax=ymax;
    model.zmax=zmax;
    model.MAPSIZE_X = MAPSIZE_X;
    model.MAPSIZE_Y = MAPSIZE_Y;
    model.X = X;
    model.Y = Y;
    model.H = H;
    model.threats =Threats;
%     PlotModel(model);
end
