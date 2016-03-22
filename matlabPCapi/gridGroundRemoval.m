function nonGroundPoints = gridGroundRemoval(lidarFrame, nGrid, cutOff)
% takes away non ground points by laying an xy-grid over the frame and
% taking away the lowest points in each grid cell based on the cutOff
% 
% in:
%   lidarFrame -    Nx3 matrix, x,y,z 
%   nGrid      -    number of grids in x- and y-direction
%   cutOff     -    percentage of a cutoff distance
%
% out:
%   nonGroundPoints - Nx3 matrix, x,y,z

nonGroundPoints = [];

p = lidarFrame;
minx = min(p(:,1));
miny = min(p(:,2));
maxx = max(p(:,1));
maxy = max(p(:,2));
gridSizeX = (maxx-minx)/nGrid;
gridSizeY = (maxy-miny)/nGrid;
zlimit = max([gridSizeX,gridSizeY])*cutOff;
for gx = 1:nGrid
    %partition into one x-grid-column
    xind = p(:,1) > minx+gridSizeX*(gx-1) & p(:,1) < minx+gridSizeX*gx;
    xp = p(xind,:);
    p(xind,:) = [];
    for gy = 1:nGrid
        %partition into all y-grid rows in that column
        yind = xp(:,2) > miny+gridSizeY*(gy-1) & xp(:,2) < miny+gridSizeY*gy;
        points = xp(yind,:);
        minz = min(points(:,3));
        zThresh = minz + zlimit;
        points(points(:,3) < zThresh,:) = [];
        if(size(points,1) > 0)
            nonGroundPoints(end+1:end+size(points,1),:) = points;
        end
    end
end