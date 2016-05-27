%% run an initial prediction with the current NN to help with tagging?
%obs: you need to train the network first to obtain the Theta values
%obs: ALWAYS SAVE YOUR CURRENT isCarMat FIRST
isCarMat = ones(Num,50);
for i=1:Num
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
        
        cluster = [pointNumber density center w l h distToEgo];
        
        %save the prediction
        %predict(Theta1, Theta2, cluster)
        cluster = (cluster-mu_nn)./sigma_nn;
        isCarMat(i,j) = predict(Theta1, Theta2, cluster);
    end
end

%% plot frame by frame, cluster by cluster to list out desired objects
figure('WindowKeyPressFcn', @keyPress)
%fClutter = @(src,evt) evalin('base','isCarMat(i,j) = 1;');
%fCar = @(src,evt) evalin('base','isCarMat(i,j) = 2;');
%hButton1 = uicontrol( 'Units', 'normalized', 'Position', [0.1 0.1 0.1 0.1], 'Style', 'pushbutton', 'Tag', 'button1', 'String', 'Clutter', 'callback', fClutter);
%hButton2 = uicontrol( 'Units', 'normalized', 'Position', [0.3 0.1 0.1 0.1], 'Style', 'pushbutton', 'Tag', 'button2', 'String', 'Car', 'callback', fCar);
%isCarMat = ones(Num,50);
i = 156;
while i <= Num
    cluster = [];
    color = [];
    for j = 1:length(clusters{i})
        pointscolor=uint8(zeros(size(clusters{i}{j},1),3));
        pointscolor(:,1)=ceil(rand(1)*50);
        pointscolor(:,2)=ceil(rand(1)*50);
        pointscolor(:,3)=ceil(rand(1)*50);
        color = [color; pointscolor];
        cluster = [cluster; clusters{i}{j}];
    end
    frameCluster = pointCloud(cluster);
    frameCluster.Color = color;
    j = 1;
    while j <= length(clusters{i})
        cluster = pointCloud(clusters{i}{j});
        pcshow(frameCluster)
        hold on
        pcshow(cluster)
        str = sprintf('Frame: %d; Cluster: %d', [i j]);
        title(str);
        testtxt = strcat('\leftarrow i: ', num2str(j));
        mu = max(clusters{i}{j});
        text(double(mu(1)), double(mu(2)), testtxt)
        %axis([150 250 50 130 60 80])
        %axis([20 200 -80 0 60 80])
        %axis([-50 50 -50 50 -5 10])
        %axis([-30 80 -30 50 -5 10])
        axis([-10 50 -50 20 -2 3])
        az = -60;
        el = 60;
        view(az, el);
        zoom(1.8)
        waitforbuttonpress;
        j = j+1;
        hold off
    end
    i = i+1;
end

%% plot frame by frame with clusternumbers assigned
%isCarMat = ones(Num,50);
figure;
i = 115;
while i <= Num
    for j = 1:length(clusters{i})
        cluster = pointCloud(clusters{i}{j});
        pcshow(cluster)
        hold on
        str = sprintf('Frame: %d', i);
        %title(str);
        text(45,-15,str)
        testtxt = strcat('\leftarrow i: ', num2str(j));
        mu = max(clusters{i}{j});
        text(double(mu(1)), double(mu(2)), testtxt) 
        axis([-10 50 -50 20 -2 3])
        az = -60;
        el = 60;
        view(az, el);
        zoom(1.8)
        j = j+1;
    end
    hold off
    waitforbuttonpress;
    i = i+1;
end

%% convert a vector of frameXclusterNumber to a class matrix
%1 for walls/clutter/noise
%isCarMat = ones(length(isCarCluster),20);
n = 114;
for i=1:n
    for j = 1:length(cars(i,:))
        if ~(cars(i,j) == 0)
            isCarMat(i,cars(i,j)) = 2;
        end
    end
    for j = 1:length(cycles(i,:))
        if ~(cycles(i,j) == 0)
            isCarMat(i,cycles(i,j)) = 3;
        end
    end
    for j = 1:length(persons(i,:))
        if ~(persons(i,j) == 0)
            isCarMat(i,persons(i,j)) = 4;
        end
    end
    for j = 1:length(groups(i,:))
        if ~(groups(i,j) == 0)
            isCarMat(i,groups(i,j)) = 5;
        end
    end
end

%% prepare ground truth cluster data for training the network
clusterObjects = [];
offset = cell(1,Num);
for i=1:Num
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
        label = isCarMat(i,j);
        clusterObjects = [clusterObjects; pointNumber density center w l h distToEgo label];
    end
end

%save data/nn_clusters_kitti_static_crossing_cyclist.mat clusterObjects
%save data/isCarMat_kitti_static_crossing_cyclist.mat isCarMat
%save data/nn_clusters_1600_1800.mat clusterObjects
%save data/isCarMat_1600_1600_1800.mat isCarMat

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
            pointscolor(:,2)=ceil(150);
            pointscolor(:,3)=ceil(150);
        elseif isCarMat(i,j) == 3
            pointscolor(:,1)=ceil(150);
            pointscolor(:,2)=ceil(255);
            pointscolor(:,3)=ceil(150);
        elseif isCarMat(i,j) == 4
            pointscolor(:,1)=ceil(150);
            pointscolor(:,2)=ceil(150);
            pointscolor(:,3)=ceil(255);
        elseif isCarMat(i,j) == 5
            pointscolor(:,1)=ceil(200);
            pointscolor(:,2)=ceil(200);
            pointscolor(:,3)=ceil(200);
        else
            pointscolor(:,1)=ceil(rand(1)*50);
            pointscolor(:,2)=ceil(rand(1)*50);
            pointscolor(:,3)=ceil(rand(1)*50);
        end
        color = [color; pointscolor];
        cluster = [cluster; clusters{i}{j}];
    end
    clusteredPC{i} = pointCloud(cluster);
    clusteredPC{i}.Color = color;
end

%% plot the clusters in their respective colors
figure
for i=1:Num
    i
    pcshow(clusteredPC{i})
    %axis([20 200 -80 0 60 80])
    axis([-10 50 -50 20 -2 3])
    zoom(2)
    pause(0.5)
end

%% cluster cars once again with a lower cutoff for Marko's corner algo
carClusters = cell(1,114);
for i=1:114
    carIndices = find(isCarMat(i,:) == 2);
    carClusters{i} = cell(1,length(carIndices));
    for j = 1:length(carIndices)        
        carClusters{i}{j} = clusters{i}{carIndices(j)};        
    end
end

%%
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
for i=1:114
    i
    for j=1:length(carClusters{i})
        pcshow(pointCloud(carClusters{i}{j}))
        hold on
        %axis([80 120 -30 -15 60 70])
        %axis([20 200 -80 0 60 80])
        axis([-10 50 -50 20 -1 2])
        zoom(2)
    end
    pause(0.1)
    hold off
end