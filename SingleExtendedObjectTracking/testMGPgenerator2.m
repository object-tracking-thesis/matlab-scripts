% Test out MGPgenerator2

close all
% Load carClusters
load carClustersCutOff.mat
load egoPosition.mat
% Get a frame,
carClusters = carClustersCutOff;
twoFrames = cell(1,2);
for N = [30 31]
    
    Ntg = carClusters{N}(:,1:3);
    Ntg(:,1:2) = Ntg(:,1:2) - repmat([egoPosition{1}(1), egoPosition{1}(2)],length(Ntg),1);
    twoFrames{1,N-29} = Ntg;
end

state30 = [-16.96, -2.28, 4.39, 4.04, 0.3, 1.74, 4.45]';

state31 = [-17.25, -2.57, 4.39, 4.07, 0.3, 1.75, 4.47]';

% [m1, m2 ,uOp] = cornerPoint(twoFrames{1})
% 
% c1 = uOp(1); c2 = uOp(2);
% n1 = uOp(3); n2 = uOp(4);
% 
% xc = (-n1*c1 + n2*c2);
% yc = (-n2*c1 -n1*c2);
% 
% xk = xc + 2.40*cos(4.04 - 0.37)
% yk = yc + 2.40*sin(4.04 - 0.37)


%

f = figure; f.Position = [100 100 1400 1000];
hold on;
p1 = plot(twoFrames{1}(:,1), twoFrames{1}(:,2),'rx');
p2 = plot([], [])%plot(twoFrames{2}(:,1), twoFrames{2}(:,2),'bx');
legend([p1 p2], 'Frame N = 30', 'Frame N = 31')
axis equal; grid on
%plot(xk, yk,'cx')
%%
nrMgps = 1;
mgpGen = MGPgenerator2(nrMgps);

clusterZ = twoFrames{1};

predictedState = state30;

corner = mgpGen.getCarCorner(clusterZ, predictedState); % This seems to be working 

[wViewed, lViewed] = mgpGen.getViewedLengths(clusterZ, predictedState); % Seems to be working as well

[orderedMgps, orderedJacobs] = mgpGen.constructMGPs(corner, predictedState, wViewed, lViewed); % This looks ok! 

assignedZ = mgpGen.assignMgps(clusterZ, predictedState); % Seems to work as well

%% Lets plot the orderedMgps and see how they stack up

figure(f)
    plot(orderedMgps(:,1), orderedMgps(:,2), 'cs')%,'MarkerFaceColor','c')
    plot(assignedZ(:,1), assignedZ(:,2), 'co')

    
    
    
    
    
    