%% cluster the frames
%clustering(pcdXYZ,max-distance,min-points-per-cluster)
clusters = cell(1,Num);
class = cell(1,Num);
for i=1:Num
    i
    tic 
    [clusters{i} class{i}] = clusterLidar(cleanedFrames{i},1,100);
    sub = cell(1,length(unique(class{i})));
    for j = 1:length(unique(class{i}))
        classes = unique(class{i});
        c = classes(j);
        ind = find(class{i} == c);
        sub{j} = clusters{i}(ind,:);
    end
    clusters{i} = sub;
    %clusters{i} = clustering(teeest,1,100);
    toc
end

%% assign different colors to all clusters found in each frame
clusteredPC = cell(1,Num);
for i=1:Num
    cluster = [];
    color = [];

    for j = 1:length(clusters{i})
        if isCarMat(i,j) == 2
            pointscolor=uint8(zeros(size(clusters{i}{j},1),3));
            pointscolor(:,1)=100;
            pointscolor(:,2)=100;
            pointscolor(:,3)=100;
            color = [color; pointscolor];
        else
            pointscolor=uint8(zeros(size(clusters{i}{j},1),3));
            pointscolor(:,1)=ceil(rand(1)*255);
            pointscolor(:,2)=ceil(rand(1)*255);
            pointscolor(:,3)=ceil(rand(1)*255);
            color = [color; pointscolor];
        end
        cluster = [cluster; clusters{i}{j}];
    end
    clusteredPC{i} = pointCloud(cluster);
    clusteredPC{i}.Color = color;
end


%%
figure
for i=1:Num
    pcshow(clusteredPC{i})
    %set(gca, 'CameraPosition', [-522.5124 -877.9707  756.3152])
    %set(gca, 'CameraViewAngle', 2.7160)
    %axis([150 250 50 130 60 80])
    axis([50 350 50 130 60 80])
    %axis([50 200 -80 0 60 80])
    zoom(2)
    pause(0.3)
end

%% samuel demo
figure
for i=1:22
    orig = pointCloud(liveFrames{i}(:,1:3));
    clean = pointCloud(cleanedFrames{i}(:,1:3));
    subplot(2,2,1)
    pcshow(orig)
    title('live view')
    zoom(2.1)
    subplot(2,2,2)
    pcshow(clean)
    title('live view, static parts removed')
    zoom(2.1)
    subplot(2,2,3)
    plot3(offset{i}(1),offset{i}(2),offset{i}(3),'.','MarkerSize',20);
    hold on
    pcshow(clusteredPC{i})
    hold off
    title('live view, clusters + ego position')
    %set(gca, 'CameraPosition', [-522.5124 -877.9707  756.3152])
    %set(gca, 'CameraViewAngle', 2.7160)
    axis([150 250 50 130 60 80])
    zoom(2.1)
    subplot(2,2,4)
    liveCloud = pointCloud(liveFrames{i}(:,1:3));
    pcshowpair(staticCloudDownsampled, liveCloud)
    title('live view vs. static map')
    zoom(2.1)
    pause(0.3)
end

%% prepare cluster data for the neural network training
clusterObjects = [];
for i=1:Num
    for j = 1:length(clusters{i})
        pointNumber = size(clusters{i}{j},1);
        maxX = max(clusters{i}{j}(:,1));
        maxY = max(clusters{i}{j}(:,2));
        minX = min(clusters{i}{j}(:,1));
        minY = min(clusters{i}{j}(:,2));
        center = [(minX+maxX)/2 (minY+maxY)/2];
        area = (maxX-minX)*(maxY-minY);
        density = pointNumber/area;
        distToEgo = sqrt(sum((center-offset{i}(1:2)').^2));
        isCar = 1;
        if isCarMat(i,j) == 2
            isCar = 2;
        end
        clusterObjects = [clusterObjects; pointNumber density center distToEgo isCar];
    end
end

save nn_clusters.mat clusterObjects

%% plot single clusters to list out the cars
for i=55:Num
    figure
    for j = 1:length(clusters{i})
        cluster = pointCloud(clusters{i}{j});
        pcshow(cluster)
        str = sprintf('Cluster: %d', [j]);
        title(str);
        axis([150 250 50 130 60 80])
        pause(2)
    end
end

%% create a matrix of which clusters in the frames are cars
isCarMat = ones(100,15);
for i=1:Num
    isCarMat(i,isCarCluster(i)) = 2;
end
isCarMat(1,1) = 1;