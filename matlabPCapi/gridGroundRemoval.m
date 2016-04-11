function lidarFrame = gridGroundRemoval(lidarFrame, nGrid, cutOff)
% takes away non ground points by laying an xy-grid over the frame and
% taking away the lowest points in each grid cell based on the cutOff
% 
% in:
%   lidarFrame -    Nx3 matrix, x,y,z 
%   nGrid      -    number of grids in x- and y-direction
%   cutOff     -    percentage of a cutoff distance
%
% out:
%   lidarFrame - Nx3 matrix, x,y,z

%settings
targetSize = [nGrid nGrid];

%prep
xx = lidarFrame(:,1);
yy = lidarFrame(:,2);
zz = lidarFrame(:,3);
minX = min(xx);
maxX = max(xx);
minY = min(yy);
maxY = max(yy);
gridSizeX = (maxX-minX)/nGrid;
gridSizeY = (maxY-minY)/nGrid;
zlimit = max([gridSizeX,gridSizeY])*cutOff;

%create the x,y cells
xxCell = (round((xx-minX)/(maxX-minX)*(targetSize(1)-1)) +1);
yyCell = (round((yy-minY)/(maxY-minY)*(targetSize(2)-1)) +1);

%create a targetSize matrix containing the min z value for each x,y
%gridcell
map = accumarray([xxCell yyCell],zz,targetSize,@(z) min(z),[]);

%build a comparison vector with each point's z coordinate in the first
%col and the respective zlimit for the point's cell in the second col
zz = [zz zeros(length(zz),1)];
for i=1:length(zz)
    zz(i,2) = map(xxCell(i),yyCell(i)) + zlimit;
end
lidarFrame(zz(:,1) < zz(:,2),:) = [];