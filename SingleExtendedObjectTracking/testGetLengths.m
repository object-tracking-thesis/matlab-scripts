% How much of each side do we see? 
              
% Load carClusters
load carClustersCutOff.mat
load egoPosition.mat
% Get a frame,
N  = 30;
carClusters = carClustersCutOff;

Ntg = carClusters{N}(:,1:3);
    
Ntg(:,1:2) = Ntg(:,1:2) - repmat([egoPosition{1}(1), egoPosition{1}(2)],length(Ntg),1);

clusterZ = Ntg;

% function [w, l]= getLengths(clusterZ, predictedState)

predictedState = [0;  % x
                  0;  % y
                  0;  % v
                  5*pi/4;  % phi (heading)
                  0;  % phiDot
                  0;  % w 
                  0]; % l


[m1,m2,uOp] = cornerPoint(clusterZ);

clusterZ = Ntg(:,1:2); 

    c1 = uOp(1); c2 = uOp(2);
    n1 = uOp(3); n2 = uOp(4);

xc = (-n1*c1 + n2*c2);
yc = (-n2*c1 -n1*c2);

% Define the two vectors spanning the L-shape, normalize them
v1 = [1, -n1/n2]; v1 = v1./norm(v1,2);
v2 = [1,  n2/n1]; v2 = v2./norm(v2,2);

% Define the heading vector
vHeading = [cos(predictedState(4)), sin(predictedState(4))];

% Find the vectors aligned with length & width of L-shape
% Orthogonal vectors have dot product ~ 0, i.e. min dot product must be the
% width vector, since it is orthogonal to the heading vector. 
if abs(dot(v1, vHeading)) < abs(dot(v2,vHeading))    
    vLength =  v2;
    vWidth  =  v1;
else
    vLength =  v1;
    vWidth  =  v2;    
end

% Project point cloud points onto each line 
top = dot(clusterZ, repmat(vWidth, length(clusterZ),1),2);
widthPoints = [top top] .* repmat(vWidth, length(clusterZ), 1);

top = dot(clusterZ, repmat(vLength, length(clusterZ),1),2);
lengthPoints = [top top] .* repmat(vLength, length(clusterZ), 1);

% Project Corners
widthCorner  = dot([xc, yc], vWidth).*vWidth;
lengthCorner = dot([xc, yc], vLength).*vLength;

% Calculate euclidian distances between corner and projected points for
% each projected line 


widthSquaredDist = (repmat(widthCorner(:,1), length(clusterZ), 1) - widthPoints(:,1)).^2 + (repmat(widthCorner(:,2), length(clusterZ), 1) - widthPoints(:,2)).^2;

[~, widthIndexMax] = max(widthSquaredDist);

lengthSquaredDist = (repmat(lengthCorner(:,1), length(clusterZ), 1) - lengthPoints(:,1)).^2 + (repmat(lengthCorner(:,2), length(clusterZ), 1) - lengthPoints(:,2)).^2;

[~, lengthIndexMax] = max(lengthSquaredDist);


W_width = norm([xc, yc] - clusterZ(widthIndexMax,:),2);

L_length = norm([xc, yc] - clusterZ(lengthIndexMax,:),2);



hold on
plot(lengthPoints(:,1), lengthPoints(:,2),'bx')
plot(lengthCorner(:,1), lengthCorner(:,2),'g*')

plot(widthPoints(:,1), widthPoints(:,2),'rx')
plot(widthCorner(:,1), widthCorner(:,2),'g*')


axis equal


%%
hold on
plot([0 v1(1)], [0 v1(2)], 'b')
plot([0 v2(1)], [0 v2(2)], 'r')
plot([0 vHeading(1)], [0 vHeading(2)], 'g')





