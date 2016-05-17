%% load data
load data/bicycleClusters.mat
n = length(bicycleClusters);
start_seq = 15;
end_seq = 35;

%% plot the clusters
figure
for i=start_seq:end_seq
    pcshow(pointCloud(bicycleClusters{i}))
    axis([0 30 0 20 -5 5])
    zoom(2)
    pause(0.5)
end

%% 