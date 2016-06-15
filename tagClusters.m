%% TAGGING: plot frame by frame, cluster by cluster to list out desired objects (good for few clusters per frame)
%isCarMat = ones(Num,50); %should be uncommented once when starting a
%tagging session
figure('WindowKeyPressFcn', @keyPress)
%check keypress.m to see what kinds of objects you can tag and how
i = 1;
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

%% TAGGING: plot frame by frame with clusternumbers assigned (better when having many clusters in a frame)
%isCarMat = ones(Num,50); %should be uncommented once when starting a
%tagging session
figure;
i = 1;
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

%% DATAPREP: prepare supervised cluster data for training the network
%once prepared, you can run objClassNN.m to train a fully connected single
%hidden layer neural network
load data/isCarMat_kitti_campus_01_186.mat %object matrix for campus example tagged by Marko&Michael
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
        clusterObjects = [clusterObjects; w l h density pointNumber/distToEgo center label];
    end
end

%save data/nn_clusters_kitti_static_crossing_cyclist.mat clusterObjects
%save data/isCarMat_kitti_static_crossing_cyclist.mat isCarMat
%save data/nn_clusters_1600_1800.mat clusterObjects
%save data/isCarMat_1600_1600_1800.mat isCarMat



%% PREDICTION: run a prediction with a previously trained network
%obs: you need to train the network first to obtain the Theta, mu_nn and sigma_nn values
%obs: ALWAYS SAVE YOUR CURRENT tagged isCarMat FIRST, because it will be
%overwritten here
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
        
        cluster = [w l h density pointNumber/distToEgo center];
        
        %save the prediction
        cluster = (cluster-mu_nn)./sigma_nn;
        isCarMat(i,j) = predict(Theta1, Theta2, cluster);
    end
end

%% CHECK TAGGING/PREDICTION: plot the objects in their respective colors
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

figure
for i=1:Num
    i
    pcshow(clusteredPC{i})
    %axis([20 200 -80 0 60 80])
    axis([-10 50 -50 20 -2 3])
    zoom(2)
    pause(0.5)
end