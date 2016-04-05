% Read all frames into memory 
path = '/Users/markocotra/Downloads/2011_09_28 2/2011_09_28_drive_0039_extract/velodyne_points/data/';
f1 = '0000000000.txt';
ff = pcd2matlab(strcat(path,f1), 0);

d = dir(path);
files = d(3:end);

nrOfFiles = length(files);
storage = cell(1,nrOfFiles);

parfor j = 1:nrOfFiles
    j
   storage{j} = pcd2matlab(strcat(path,files(j).name),0);
   storage{j} = storage{j}(:,1:3); % Keep only relevant values
   storage{j} = pointCloud(storage{j}); % Make point cloud object 
end

%% 

ptCloud = storage;
objStorage = cell(1,length(ptCloud));
for k = 1:3%length(ptCloud)
    h = warndlg('Press "ok" to activate ginput', 'Activate ginput');
    h.Position = [1550,100, 200, 100];
    s = sprintf('k = %d',k);
    f = figure('name',s,'Position', [0 0 1600 1200]);
    plot(ptCloud{k}.Location(:,1),...
        ptCloud{k}.Location(:,2),'.');
    text(15,15,'PRESS ENTER TO QUIT GINPUT');
    ax = gca;
    ax.Position = [0 0 1 1];
    axis([-50 50 -50 50]);
    waitfor(h)
    h = 1;
    while h
        try
            [x, y] = ginput(1);
            r = 0.5;            
            phi = 0:0.01:2*pi;
            X = x + r*cos(phi);
            Y = y + r*sin(phi);
            hold on
            plot(X,Y,'r')
            objStorage{k} = [objStorage{k} {X;Y}];
        catch e
            h = 0;
        end
    end
    close(f)
end

objStorage = objStorage(~cellfun(@isempty, objStorage));
%% Filter out the markings 
[~,c] = size(objStorage);
keepPoints = cell(1,c);

for k = 1:c    
    points = ptCloud{k}.Location;
    for h = 1:size((objStorage{k}),2)
        
        xmin = min(objStorage{k}{1,h});
        xmax = max(objStorage{k}{1,h});
        ymin = min(objStorage{k}{2,h});
        ymax = max(objStorage{k}{2,h});
        
        tPoints = points(points(:,1) > xmin,:,:);
        tPoints = tPoints(tPoints(:,1) < xmax,:,:);
        tPoints = tPoints(tPoints(:,2) > ymin,:,:);
        tPoints = tPoints(tPoints(:,2) < ymax,:,:);
       
        keepPoints{k} = [keepPoints{k} {pointCloud(tPoints)}];
    end
end
%% Test to plot 








