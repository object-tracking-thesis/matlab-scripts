%% test points
points = [1 1; 0.9 0.9; 1.4 1; 3 3; 2.9 2.9; 1.1 1.1];

%% real points
%use ground-removed points here
points = lidarData{8};

%% find clusters
%clustering(pcdXYZ,max-distance,min-points-per-cluster)
clusters = clustering(points,0.5,50);

%% plot all clusters
figure;
wd = 50;

for i = 1:length(clusters)
    h = lidarPlot(clusters{i});
    axis(wd*[-1 1 -1 1 -1/(0.5*wd) 1])
    grid off; box off;  
    color = rand(1,3);
    set(h,'MarkerEdgeColor', color,'MarkerFaceColor', color);
    hold on
end