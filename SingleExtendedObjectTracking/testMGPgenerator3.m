% test MGPgenerator3 

load carClustersCutOff.mat
load egoPosition.mat
% Get a frame,
carClusters = carClustersCutOff;

Ntg = carClusters{30}(:,1:3);
Ntg(:,1:2) = Ntg(:,1:2) - repmat([egoPosition{1}(1), egoPosition{1}(2)],length(Ntg),1);
  

clusterZ = Ntg;
%%

state30 = [-16.96, -2.28, 4.39, 4.0, 0.3, 1.6, 4.45]';

    drawMyRide(state30,'b');
        axis equal; hold on 
    plot(clusterZ(:,1), clusterZ(:,2),'kx')

mgpGen3 = MGPgenerator3(0);

corner = mgpGen3.getCarCorner(clusterZ, state30);

[wViewed, lViewed] = mgpGen3.getViewedLengths(clusterZ, state30);

assignedZ = mgpGen3.assignMgps(clusterZ, state30);

    plot(assignedZ(:,1), assignedZ(:,2), 'rx')

mgpHandles = mgpGen3.constructMGPs(corner, state30, wViewed, lViewed);

MGP1 = mgpHandles{1}(state30);
MGP2 = mgpHandles{2}(state30);
MGP3 = mgpHandles{3}(state30);

plot(MGP1(1), MGP1(2), 'oc')
plot(MGP2(1), MGP2(2), 'oc')
plot(MGP3(1), MGP3(2), 'oc')

