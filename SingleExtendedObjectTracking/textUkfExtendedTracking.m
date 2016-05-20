% textUkfExtendedTracking.m

load carClustersCutOff.mat
load egoPosition.mat
% Get a frame,
carClusters = carClustersCutOff;

Ntg = carClusters{30}(:,1:3);
Ntg(:,1:2) = Ntg(:,1:2) - repmat([egoPosition{1}(1), egoPosition{1}(2)],length(Ntg),1);  

clusterZ30 = Ntg;

Ntg = carClusters{31}(:,1:3);
Ntg(:,1:2) = Ntg(:,1:2) - repmat([egoPosition{1}(1), egoPosition{1}(2)],length(Ntg),1);  

clusterZ31 = Ntg;

%% Test out UKF against state30 --> state31, 
%  where state30 is x0 and the measurements at 31 are from the next state

% st = [x_k, y_k, v_k, phi_k, phiDot w_k, l_k]';
T = 0.1; % sample time 
f = @(st) [st(1)+T*st(3)*cos(st(4));...
           st(2)+T*st(3)*sin(st(4));...
           st(3);...
           st(4) + T*st(5);...
           st(5);
           st(6);
           st(7)];

% Define motion noise
velCov = 0.05^2;
phiDotCov = 0.05^2;
% Define width and length noise 
wCov = 0.01^2; % std.dev is 3 cm
lCov = 0.02^2; % std.dev is 5 cm
R = 0.3^2;
subQ = diag([velCov, phiDotCov, wCov, lCov]);
gamma = [0 0 1 0 0 0 0 ;
         0 0 0 0 1 0 0 ;
         0 0 0 0 0 1 0 ;
         0 0 0 0 0 0 1]';
 
% MOTION COVARIANCE MATRIX 
Q = T*gamma*subQ*gamma';
 
% Define base state & base covariance (frame N = 30 in carClustersCutOff) 
st30 = [-16.96, -2.28, 4.39, 4, 0, 1.74, 4.45]';

st31 = [-17.25, -2.57, 4.41, 4.05, 0.31, 1.75, 4.47]';

P0 = cov([st30';
          st31']);
P0 = P;      
%

%plot(clusterZ31(:,1), clusterZ31(:,2),'kx'); hold on
%drawMyRide(st31,'c'); axis equal

% Initialize UKF and use it to predict 
nMGPS = 7;
nObsSt = 2;
nSt = 7;
ukf = UKF(Q,R, nMGPS, nObsSt, nSt, st30, P0);
ukf.predictMoments(f);

predictedState = ukf.predSt;

% Initialize mgpGen3 and use it to create MGPs
N = 2;
mgpGen3 = MGPgenerator3(N);

               corner = mgpGen3.getCarCorner(clusterZ31, predictedState);
   [wViewed, lViewed] = mgpGen3.getViewedLengths(clusterZ31, predictedState);
            assignedZ = mgpGen3.assignMgps(clusterZ31, predictedState);

         [mgpHandles] = mgpGen3.constructMGPs(corner, predictedState, wViewed, lViewed);



% Update UKF
assignedZo = reshape(assignedZ', 14,1);
%ukf.setNrMGPS(length(assignedZ));
ukf.updateMoments(mgpHandles, assignedZo);
ukf.upSt

% Let's see how we did 
fig = figure; fig.Position = [100 100 1200 800];
plot(clusterZ31(:,1), clusterZ31(:,2),'rx');hold on; axis equal
%drawMyRide(st30,'k'); 
drawMyRide(st31,'k','--')
drawMyRide(predictedState, 'c')
drawMyRide(ukf.upSt,'b')

plot(assignedZ(:,1), assignedZ(:,2),'g*','MarkerSize',10)

for j = 1:length(mgpHandles)
   A = mgpHandles{j}(predictedState); 
   plot(A(1), A(2),'ok')
end






