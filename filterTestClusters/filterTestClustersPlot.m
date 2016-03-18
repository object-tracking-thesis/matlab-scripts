%basic plot of the filterTestClusters

%the straight line in direction north-west is the car going
%away from us, all the rest is clutter that sometimes appears or disappears
%(it's usually walls, sometimes they land above and sometimes below the
%current clustering threshold (100points) and therefore seem to "appear"
%and "disappear")
%they also slightly change their position from frame to frame due to the
%fact that I currently just compute the mean of each cluster which varies
%slightly from frame to frame

%% load the data
%the clusters are cells of xy-column matrices, so if there are 4 clusters
%in a frame, there will be a 4x2 matrix in that cell

%filterTestClustersXY1 -> following behind another car into the crossing, 
%the other car then proceeds straight ahead while we wait at the crossing
load('filterTestClustersXY1.mat', 'filterTestClusters')

%filterTestClustersXY2 -> following behind another car into the crossing,
%the other car then turns left while we wait at the crossing
%load('filterTestClustersXY2.mat', 'filterTestClusters')

%% plot to check that everything went alright
figure
for i=1:length(filterTestClusters)
    scatter(filterTestClusters{i}(:,1), filterTestClusters{i}(:,2),'x');
    hold on
    pause(0.5)
end
hold off

