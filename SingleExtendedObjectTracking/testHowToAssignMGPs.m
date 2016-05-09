% Test out how to assign MGPS

% Load carClusters
load carClustersCutOff.mat
load egoPosition.mat
% Get a frame,
N  = 30;
carClusters = carClustersCutOff;

Ntg = carClusters{N}(:,1:3);
predictedState = [0 0 0 5*pi/4 0 0 0]';
    
Ntg(:,1:2) = Ntg(:,1:2) - repmat([egoPosition{1}(1), egoPosition{1}(2)],length(Ntg),1);

[m1,m2,uOp, filtNtg] = cornerPoint(Ntg); 

    c1 = uOp(1); c2 = uOp(2);
    n1 = uOp(3); n2 = uOp(4);

xc = (-n1*c1 + n2*c2);
yc = (-n2*c1 -n1*c2);

    
M = mean(Ntg(:,1:2));

% Formulate all 4 possible L-vectors, to find the ones that are along the
% data

v1 = [1, -n1/n2]; v1 = v1./norm(v1,2);
v2 = [1,  n2/n1]; v2 = v2./norm(v2,2);
v3 = -1.*v1;
v4 = -1.*v2;

p1 = [xc yc] + v1;
p2 = [xc yc] + v2;
p3 = [xc yc] + v3;
p4 = [xc yc] + v4;

d1 = (p1(1) - M(1))^2 + (p1(2) - M(2))^2;
d2 = (p2(1) - M(1))^2 + (p2(2) - M(2))^2;
d3 = (p3(1) - M(1))^2 + (p3(2) - M(2))^2;
d4 = (p4(1) - M(1))^2 + (p4(2) - M(2))^2;

% The massVec* are vectors which are in the direction of the point cloud
% points, from [xc yc]

if d1 < d3
    massVec1 = v1;
else
    massVec1 = v3;
end

if d2 < d4
    massVec2 = v2;
else
    massVec2 = v4;
end


% Define the heading vector
vHeading = [cos(predictedState(4)), sin(predictedState(4))];

% Find the vectors aligned with length & width of L-shape
% Orthogonal vectors have dot product ~ 0, i.e. min dot product must be the
% width vector, since it is orthogonal to the heading vector.
if abs(dot(massVec1, vHeading)) < abs(dot(massVec2,vHeading))
    vLength =  massVec2;
    vWidth  =  massVec1;
else
    vLength =  massVec1;
    vWidth  =  massVec2;
end

N = 1; % Test assumption

% Get the viewed lengths of each side 

% [wViewed, lViewed] = getViewedLengths(~, clusterZ, predictedState)

wViewed = 1.6016;
lViewed = 4.2166;

storageCP = zeros(1,2); % Storage for the cornerPoint 
storageCPmeas = zeros(1,2); % Storage for the assiociated measurementPoint to cornerPoint
storageCP = [xc, yc];

storageW = [];
storageWmeas = [];
storageL = [];
storageLmeas = [];

if wViewed > 0.5         
    storageW = zeros(N+1,2); % Storage for the practical Width MGPs
    storageWmeas = zeros(N+1,2);
    
    for j = 1:N+1
        storageW(j,:) = [xc, yc] + wViewed/(N+1)*j*vWidth;
    end
end

if lViewed > 0.5
    storageL = zeros(N+1,2); % Storage for the practical Length MGPs
    storageLmeas = zeros(N+1,2);    
    
    for j = 1:N+1
        storageL(j,:) = [xc, yc] + lViewed/(N+1)*j*vLength;
    end
end

allStorage = [storageL; 
              storageCP; 
              storageW];                    
          
IDX = knnsearch(Ntg(:,1:2), allStorage);

chosenMeasurements = Ntg(IDX,1:2);













