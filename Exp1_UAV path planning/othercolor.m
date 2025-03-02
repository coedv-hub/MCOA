function c = othercolor(n,m)

types = who('-file','colorData.mat');

% if no colormap is choosen then display available colormaps
if nargin < 1
    c = types;
else
    % default number of points
    if nargin < 2, m = size(get(gcf,'colormap'),1); end

    % allows numerical indexing
    if isnumeric(n), n = char(types(n)); end
        
    % load color data
    data = load('colorData.mat',n);
    if isempty(fieldnames(data))
        c = [];
    else
        c = interp1(linspace(0,1,size(data.(n),1)),data.(n),linspace(0,1,m),'cubic');
        c(c<0) = 0;
        c(c>1) = 1;
    end
end