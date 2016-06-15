%% cluster the frames
clusters = cell(1,Num);
class = cell(1,Num);
for i=1:Num
    i
    tic 
    [clusters{i} class{i}] = clusterLidar(cleanedFrames{i},0.7,50);
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
    zoom(1.4)
    axis([-20 50 -40 20 -4 4])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    pause(0.3)
end