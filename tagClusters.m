%% run an initial prediction with the current NN to help with tagging?
%obs: you need to train the network first to obtain the Theta values
isCarMat = ones(Num,20);
for i=1:150
    for j = 1:length(clusters{i})
        pointNumber = size(clusters{i}{j},1);
        maxX = max(clusters{i}{j}(:,1));
        maxY = max(clusters{i}{j}(:,2));
        maxZ = max(clusters{i}{j}(:,3));
        minX = min(clusters{i}{j}(:,1));
        minY = min(clusters{i}{j}(:,2));
        minZ = min(clusters{i}{j}(:,3));
        center = [(minX+maxX)/2 (minY+maxY)/2 (minZ+maxZ)/2];
        w = (maxX-minX);
        l = (maxY-minY);
        h = (maxZ-minZ);
        volume = (maxX-minX)*(maxY-minY)*(maxZ-minZ);
        density = pointNumber/volume;
        distToEgo = sqrt(sum((center-offset{i}(1:3)').^2));
        
        cluster = [pointNumber density center w l h distToEgo];
        
        %save the prediction
        %predict(Theta1, Theta2, cluster)
        cluster = (cluster-mu)./sigma;
        isCarMat(i,j) = predict(Theta1, Theta2, cluster);
    end
end

%% plot frame by frame, cluster by cluster to list out desired objects
figure('WindowKeyPressFcn', @keyPress)
%fClutter = @(src,evt) evalin('base','isCarMat(i,j) = 1;');
%fCar = @(src,evt) evalin('base','isCarMat(i,j) = 2;');
%hButton1 = uicontrol( 'Units', 'normalized', 'Position', [0.1 0.1 0.1 0.1], 'Style', 'pushbutton', 'Tag', 'button1', 'String', 'Clutter', 'callback', fClutter);
%hButton2 = uicontrol( 'Units', 'normalized', 'Position', [0.3 0.1 0.1 0.1], 'Style', 'pushbutton', 'Tag', 'button2', 'String', 'Car', 'callback', fCar);
%isCarMat = ones(Num,20);
i = 89;
while i <= Num
    j = 1;
    while j <= length(clusters{i})
        cluster = pointCloud(clusters{i}{j});
        pcshow(cluster)
        str = sprintf('Frame: %d; Cluster: %d', [i j]);
        title(str);
        %axis([150 250 50 130 60 80])
        axis([20 200 -80 0 60 80])
        %axis([-50 50 -50 50 -5 10])
        zoom(1.5)
        waitforbuttonpress;
        j = j+1;
    end
    i = i+1;
end

%% convert a vector of frameXclusterNumber to a class matrix
%1 for walls/clutter/noise
isCarMat = ones(length(isCarCluster),20);
for i=1:size(isCarMat,1)
    %99 if this frame does not contain any desired object
    if isCarCluster(i) == 99
        continue
    end
    %2 for car
    isCarMat(i,isCarCluster(i)) = 2;
end

%% prepare ground truth cluster data for training the network
clusterObjects = [];
offset = cell(1,Num);
for i=1:200
    offset{i} = [0 0 0]';
    for j = 1:length(clusters{i})
        pointNumber = size(clusters{i}{j},1);
        maxX = max(clusters{i}{j}(:,1));
        maxY = max(clusters{i}{j}(:,2));
        maxZ = max(clusters{i}{j}(:,3));
        minX = min(clusters{i}{j}(:,1));
        minY = min(clusters{i}{j}(:,2));
        minZ = min(clusters{i}{j}(:,3));
        center = [(minX+maxX)/2 (minY+maxY)/2 (minZ+maxZ)/2];      
        w = (maxX-minX);
        l = (maxY-minY);
        h = (maxZ-minZ);
        volume = (maxX-minX)*(maxY-minY)*(maxZ-minZ);
        density = pointNumber/volume;
        distToEgo = sqrt(sum((center-offset{i}(1:3)').^2));
        isCar = 1;
        if isCarMat(i,j) == 2
            isCar = 2;
        end
        clusterObjects = [clusterObjects; pointNumber density center w l h distToEgo isCar];
    end
end

%save data/nn_clusters_kitti_static_crossing.mat clusterObjects
%save data/isCarMat_kitti_static_crossing.mat isCarMat
save data/nn_clusters_1600_1800.mat clusterObjects
save data/isCarMat_1600_1600_1800.mat isCarMat

%% assign different colors to all clusters found in each frame
clusteredPC = cell(1,Num);
for i=1:Num
    cluster = [];
    color = [];
    for j = 1:length(clusters{i})
        pointscolor=uint8(zeros(size(clusters{i}{j},1),3));
        %assign the same color for all cars
        if isCarMat(i,j) == 2
            pointscolor(:,1)=ceil(255);
            pointscolor(:,2)=ceil(200);
            pointscolor(:,3)=ceil(200);
        else
            pointscolor(:,1)=ceil(rand(1)*150);
            pointscolor(:,2)=ceil(rand(1)*150);
            pointscolor(:,3)=ceil(rand(1)*150);
        end
        color = [color; pointscolor];
        cluster = [cluster; clusters{i}{j}];
    end
    clusteredPC{i} = pointCloud(cluster);
    clusteredPC{i}.Color = color;
end

%% plot the clusters in their respective colors
figure
for i=1:150
    i
    pcshow(clusteredPC{i})
    axis([20 200 -80 0 60 80])
    zoom(2)
    pause(0.1)
end

%% cluster cars once again with a lower cutoff for Marko's corner algo
carClusters = cell(1,150);
for i=1:150
    for j = 1:length(clusters{i})
        if isCarMat(i,j) == 2
            carClusters{i} = clusters{i}{j};
        end
    end
end
carClustersCutOff = cell(1,150);
classCutOff = cell(1,150);
for i=1:150
    i
    tic 
    [carClustersCutOff{i}, classCutOff{i}] = clusterLidar(carClusters{i},0.2,50);
    sub = cell(1,length(unique(classCutOff{i})));
    for j = 1:length(unique(classCutOff{i}))
        classes = unique(classCutOff{i});
        c = classes(j);
        ind = find(classCutOff{i} == c);
        sub{j} = carClustersCutOff{i}(ind,:);
    end
    carClustersCutOff{i} = sub{1};
    toc
end
%%
figure;
for i=1:150
    i
    pcshow(pointCloud(carClustersCutOff{i}{1}))
    %axis([80 120 -30 -15 60 70])
    axis([20 200 -80 0 60 80])
    zoom(2)
    pause(0.1)
end