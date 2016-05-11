% Filter data to make it better, using cP()

% Load carClusters
load carClustersCutOff.mat
load egoPosition.mat
% Get a frame,
N  = 30;
carClusters = carClustersCutOff;

Ntg = carClusters{N}(:,1:3);
    
Ntg(:,1:2) = Ntg(:,1:2) - repmat([egoPosition{1}(1), egoPosition{1}(2)],length(Ntg),1);

[m1,m2,uOp, filtNtg] = cornerPoint(Ntg); % 0.3 0.5

    c1 = uOp(1); c2 = uOp(2);
    n1 = uOp(3); n2 = uOp(4);

xc = (-n1*c1 + n2*c2);
yc = (-n2*c1 -n1*c2);
    
    
f = figure;
f.Position = [100 100 1200 800];

plot(filtNtg(:,1), filtNtg(:,2),'kx'); hold on 
plot(xc, yc,'c*')
%plot(filtNtg(23,1), filtNtg(23,2),'rx')
%plot(filtNtg(488,1), filtNtg(488,2),'gx')
M1 = mean(filtNtg);
M2 = median(filtNtg);

%plot(M1(1), M1(2),'rx')
%plot(M2(1), M2(2),'gx')

% ------
% This code requires the scripts: testHowToAssignMGPs.m & testGetLengths.m
% ------

p1 = plot(Ntg(widthIndexMax,1), Ntg(widthIndexMax,2),'gs','MarkerSize',10);
p2 = plot(Ntg(lengthIndexMax,1), Ntg(lengthIndexMax,2),'gs','MarkerSize',10);

axis equal

p3 = plot(storageW(:,1), storageW(:,2),'om');
p4 = plot(storageL(:,1), storageL(:,2),'oc');

p5 = plot(chosenMeasurements(:,1), chosenMeasurements(:,2), 'rx');

legend([p1 p2 p3 p4 p5],...
       'Maximum spanning width point',...
       'Maximum spanning length point',...
       'practical MGPs for width',...
       'practical MGPs for length',...
       'Chosen measurements')



predictedState = [0 0 0 5*pi/4 0];

getCarCorner(Ntg, predictedState)
