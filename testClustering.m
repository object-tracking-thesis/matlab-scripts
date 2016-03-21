%% cluster the frames
%clustering(pcdXYZ,max-distance,min-points-per-cluster)
clusters = cell(1,Num);
for i=1:Num
    i
    tic 
    clusters{i} = clustering(cleanedFrames{i},1,100);
    %clusters{i} = clustering(teeest,1,100);
    toc
end

%% assign different colors to all clusters found in each frame
clusteredPC = cell(1,Num);
for i=1:Num
    cluster = [];
    color = [];

    for j = 1:length(clusters{i})
        pointscolor=uint8(zeros(size(clusters{i}{j},1),3));
        pointscolor(:,1)=ceil(rand(1)*255);
        pointscolor(:,2)=ceil(rand(1)*255);
        pointscolor(:,3)=ceil(rand(1)*255);
        color = [color; pointscolor];
        cluster = [cluster; clusters{i}{j}];
    end
    clusteredPC{i} = pointCloud(cluster);
    clusteredPC{i}.Color = color;
end


%%
figure
for i=40:90
    pcshow(clusteredPC{i})
    %set(gca, 'CameraPosition', [-522.5124 -877.9707  756.3152])
    %set(gca, 'CameraViewAngle', 2.7160)
    axis([150 250 50 130 60 80])
    zoom(2)
    pause(0.3)
end

%% samuel demo
figure
for i=70:120
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