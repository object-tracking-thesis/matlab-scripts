%% cluster the frames
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
    fprintf('number of clusters: %6.2d\n', length(clusters{i}));
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


%% plot the clusters in their respective colors
figure
for i=1:Num
    i
    pcshow(clusteredPC{i})
    hold on
    %plot3(egoPosition{i}(1),egoPosition{i}(2),egoPosition{i}(3),'x')
    %axis([20 200 -80 0 60 80])
    zoom(2)
    pause(0.3)
    hold off
end