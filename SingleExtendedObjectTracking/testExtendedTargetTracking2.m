% test Extended Target Tracking using an Extended KF with the new
% MGPgenerator2 class. Investigated for timesteps N = 30 to 31, with the
% carClusters mat-file. 

% Most of the stuff is directly copied from textExtendedTargetTracking.m
%%
close all
 
% Define nonlinear motion model 
%       1    2    3    4      5      6    7
% st = [x_k, y_k, v_k, phi_k, phiDot w_k, l_k]';
T = 0.1; % sample time 
f = @(st) [st(1)+T*st(3)*cos(st(4));...
           st(2)+T*st(3)*sin(st(4));...
           st(3);...
           st(4) + T*st(5);...
           st(5);
           st(6);
           st(7)];

% Define Jacobian for nonlinear model
J = @(st) [[ 1, 0, T*cos(st(4)),  -T*st(3)*sin(st(4)), 0, 0, 0];
           [ 0, 1, T*sin(st(4)),   T*st(3)*cos(st(4)), 0, 0, 0];
           [ 0, 0,            1,                    0, 0, 0, 0];
           [ 0, 0,            0,                    1, T, 0, 0];
           [ 0, 0,            0,                    0, 1, 0, 0];
           [ 0, 0,            0,                    0, 0, 1, 0];
           [ 0, 0,            0,                    0, 0, 0, 1]];

% Define motion noise
velCov = 1^2;
phiDotCov = 0.1^2;
% Define width and length noise 
wCov = 0.03^2; % std.dev is 3 cm
lCov = 0.05^2; % std.dev is 5 cm
 
subQ = diag([velCov, phiDotCov, wCov, lCov]);
gamma = [0 0 1 0 0 0 0 ;
         0 0 0 0 1 0 0 ;
         0 0 0 0 0 1 0 ;
         0 0 0 0 0 0 1]';
 
% MOTION COVARIANCE MATRIX 
Q = T*gamma*subQ*gamma';
 
% Define base state & base covariance (frame N = 30 in carClustersCutOff) 
st30 = [-16.96, -2.28, 4.39, 4.04, 0.29, 1.74, 4.45]';
st31 = [-17.25, -2.57, 4.41, 4.07, 0.31, 1.75, 4.47]';

P0 = cov([st30';
          st31']);


%% Load carClusters
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

% twoFrames have the clusterZ data for N = 30 to N = 31

%% Initiate a MGPgenerator2
nrMgps = 1;
mgpGen2 = MGPgenerator2(nrMgps);
%% 

%% Write Update function    
nrIter = 1;
 
    % Prediction
    predictedState = f(st30); 
    predictedCov = J(st30)*P0*J(st30)' + Q;        
   
    % Update
    corner = mgpGen2.getCarCorner(twoFrames{2}, predictedState(:,1));
    
    [wViewed, lViewed] = mgpGen2.getViewedLengths(twoFrames{2}, predictedState);
    [orderedMgps, orderedJacobs] = mgpGen2.constructMGPs(corner, predictedState, wViewed, lViewed);
    assignedZ = mgpGen2.assignMgps(twoFrames{2}, predictedState);
    
    
    
    % Stack them vertically
    len = 2*size(assignedZ,1); % How many rows do we have
    vMgps = reshape(orderedMgps.',len,1);
    vZ = reshape(assignedZ.',len,1);
    
    vJacobs = reshape(permute(orderedJacobs, [1,3,2]),[],size(orderedJacobs,2));
    
    S = vJacobs*predictedCov*vJacobs' + (0.1^2)*eye(len);
    K = predictedCov*vJacobs'*inv(S);
    
    updatedState = predictedState + K*(vZ - vMgps);
    updatedCov = predictedCov - K*S*K';

%%

plot(orderedMgps(:,1), orderedMgps(:,2),'kx')
hold on
plot(assignedZ(:,1), assignedZ(:,2), 'bo')
axis equal
hold on
for k = 1:3+nrMgps*2
   plot([orderedMgps(k,1) assignedZ(k,1)], [orderedMgps(k,2) assignedZ(k,2)], '-r')    
end


%% Let's see how the filter behaves! 

f = figure; f.Position = [100 100 1200 800];
    hold on; grid on; op = 0.6;
    drawMyRide(st30,op.*ones(1,3), ':')
    drawMyRide(predictedState,'c')
    drawMyRide(updatedState,'b')
    drawMyRide(st31,op.*ones(1,3), ':')
    
    axis equal




















        