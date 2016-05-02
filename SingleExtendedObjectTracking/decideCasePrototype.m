%
% prototype for decideCase function, needs to decide which of the following
% cases are in effect for the pointCloud: 
% 
% y CASE 1    y  CASE 2   y CASE 3
% |      o    |    o      |    
% |     /     |   /       | o   
% | o  /      |  /        |  \  
% |  \o       | o         |   \  
% |           |  \o       |    \o
% |________ x |________ x |________ x
%
% Currently, case 3 is omitted, since we perhaps don't even need it 

% Load data for one frame and pre-process it 

close all; clc
% Load carClusters
load carClusters.mat
load egoPosition.mat
% Get a frame,
N  = 84; % 40 is good
Ntg = carClusters{N}(:,1:3);

%

f = figure; f.Position = [100 100 1400 600];
    subplot(1,2,1)
    hold on; grid on; op = 0.8;
    p1 = plot(Ntg(:,1), Ntg(:,2),'x','Color',op.*[1 1 1]);
    %axis([-90 90 -90 90])
    axis equal

[m1,m2,uOp, filtNtg] = cornerPoint(Ntg, 0.1, 0);

figure(f)
    subplot(1,2,1)
    p2 = plot(filtNtg(:,1), filtNtg(:,2),'x','Color','k');
    legend([p1 p2], 'Raw Points', 'Filtered Points')
    
    subplot(1,2,2)    
    histogram(Ntg(:,3))
    title('Distribution of height measurements')
    
%
% Make copy of filtNtg
copyNtg = filtNtg; 
% Find max length vector 
%
% ====  TODO 
% This should be changed. We don't actually want max length, what we want 
% is the vector that is aligned with the velocity.
%
c1 = uOp(1); c2 = uOp(2);
n1 = uOp(3); n2 = uOp(4);

xc = (-n1*c1 + n2*c2); % Corner point
yc = (-n2*c1 - n1*c2); % Corner point 
    subplot(1,2,1)
    plot(xc, yc, 'mx')

CP = [xc, yc];

V1 = [1 -n1/n2]./norm([1 -n1/n2],2); % unitvector For m1 measure
V2 = [1  n2/n1]./norm([1  n2/n1],2); % unitvector For m2 measure 
V = [V1; V2];

[my, i] = max([m1 m2]);
V = V(i,:); % This is our max length vector

% Find rotation angle from origin
Yline = [0 1];

phi = acos(dot(V,Yline)); % Finds angle between (0,1) vector and V

%phi = mod(pi,phi);

R = [cos(phi) -sin(phi); sin(phi) cos(phi)];

% Put L-shape in origin, compensate rotation 

copyNtg = copyNtg - repmat(mean(filtNtg), length(copyNtg), 1);
CP = CP - repmat(mean(filtNtg), size(CP,1), 1);
%
copyNtg = copyNtg*R';
CP = CP*R';

figure 
plot(copyNtg(:,1), copyNtg(:,2),'x'); axis equal; hold on; grid on

plot(CP(1), CP(2), 'mx')

if CP(1)*CP(2) > 0
   title('Case 2: Back Left')
else
   title('Case 1: Back Right')
end

%% Try another approach
% Locate corner point and from that determine the vectors from (xc,yc)
% along the sides. In order to determine case, get velocity vector and look
% whether the majority is either right or left of velocity vector. 

N  = 40; % 40 is good
Ntg = carClusters{N}(:,1:3);

f = figure; f.Position = [100 100 1400 600];
    subplot(1,2,1)
    hold on; grid on; op = 0.8;
    p1 = plot(Ntg(:,1), Ntg(:,2),'x','Color',op.*[1 1 1]);
    %axis([-90 90 -90 90])
    axis equal

[m1,m2,uOp, filtNtg] = cornerPoint(Ntg, 0.1, 0);

figure(f)
    subplot(1,2,1)
    p2 = plot(filtNtg(:,1), filtNtg(:,2),'x','Color','k');
    legend([p1 p2], 'Raw Points', 'Filtered Points')
    
    subplot(1,2,2)    
    histogram(Ntg(:,3))
    title('Distribution of height measurements')



c1 = uOp(1); c2 = uOp(2);
n1 = uOp(3); n2 = uOp(4);

xc = (-n1*c1 + n2*c2); % Corner point
yc = (-n2*c1 - n1*c2); % Corner point 
    
angle = ones(1,4)*atan2(n2,n1) + [0 1 2 3]*pi/2;

m = mean(Ntg(:,1:2));
mx = m(1); my = m(2);

V = zeros(4,2);

for j = 1:4
   V(j,:) = [cos(angle(j)), sin(angle(j))] + [xc ,yc];
end

D = zeros(4,1);

for j = 1:4
   D(j,:) = sqrt((V(j,1) - mx)^2 + (V(j,2) - my)^2);
end

[~, i] = min(D);
V1 = V(i,:); % One of the vectors containing points along xc,yc
D(i) = [];
[~, i] = min(D); 
V2 = V(i,:); % The other vector containing points along xc, yc

% Place points in origin again
V1origin = V1 - [xc, yc];
V1origin = V1origin./norm(V1origin,2); % Normalize 
V2origin = V2 - [xc, yc];
V2origin = V2origin./norm(V2origin,2); % Normalize 


m1vector = [1 -n1/n2]./norm([1 -n1/n2],2); 
m2vector = [1 n2/n1]./norm([1 n2/n1],2);

if dot(m1vector, V1origin) > 0.9 % Perpendicular 
    side1 = [V1, m2]; % side1 is for vector 1
    side2 = [V2, m1]; % side2 is for vector 2
else
    side1 = [V1, m1];
    side2 = [V2, m2];
end
    
    

figure(f)
    subplot(1,2,1)
        plot(xc,yc,'ms')
        plot(V1(1), V1(2),'rs','MarkerSize',20,'MarkerFaceColor','r')
        plot(V2(1), V2(2),'bs','MarkerSize',20,'MarkerFaceColor','b')
        
        plot(V(4,1) , V(4,2),'gx')

%% Let's say we know the velocity angle, for N = 84 it's angle(2)
phi = angle(2);

% construct velocity vector 






