clear all
load carClustersCutOff.mat
load egoPosition.mat
% Get a frame,
carClusters = carClustersCutOff;

clusterZ = cell(1,30);

for j = 1:30

    Ntg = carClusters{30+j}(:,1:3);
    Ntg(:,1:2) = Ntg(:,1:2) - repmat([egoPosition{1}(1), egoPosition{1}(2)],length(Ntg),1);  

    clusterZ{j} = Ntg;
end

%% 
for j = 1:30    
   plot(clusterZ{j}(:,1),clusterZ{j}(:,2),'x'); hold on; axis equal
   pause(0.2)    
end
%%
% st = [x_k, y_k, v_k, phi_k, phiDot w_k, l_k]';
T = 0.1; % sample time 
f = @(st) [st(1)+T*st(3)*cos(st(4));...
           st(2)+T*st(3)*sin(st(4));...
           st(3);...
           st(4) + T*st(5);...
           st(5);
           st(6);
           st(7)];

run('mcSimOfP0.m')       
       
% Define motion noise
velCov = 0.5^2;    % 0.05 
phiDotCov = 0.05^2; % 0.05
% Define width and length noise 
wCov = 0.01^2; % std.dev is 1 cm
lCov = 0.02^2; % std.dev is 2 cm
R = 0.3^2;          % 0.3
subQ = diag([velCov, phiDotCov, wCov, lCov]);
gamma = [0 0 1 0 0 0 0 ;
         0 0 0 0 1 0 0 ;
         0 0 0 0 0 1 0 ;
         0 0 0 0 0 0 1]';
 
% MOTION COVARIANCE MATRIX 
Q = T*gamma*subQ*gamma';

x0 = [-16.96, -2.28, 4.39, 4, 0, 1.74, 4.45]';
P0 = P;  

%% Initialize UKF and use it to predict

for k = 1:30
    if k == 1
        nMGPS = 7;
        nObsSt = 2;
        nSt = 7;

        ukf = UKF(Q,R, nMGPS, nObsSt, nSt, x0, P0);

        ukf.predictMoments(f);
        predictedState = ukf.predSt;
    else
        ukf.predictMoments(f);
        predictedState = ukf.predSt;    
    end

    % Initialize mgpGen3 and use it to create MGPs
    if k == 1
        N = 2;
        mgpGen3 = MGPgenerator3(N);
    end     
        corner = mgpGen3.getCarCorner(clusterZ{k}, predictedState);
        [wViewed, lViewed] = mgpGen3.getViewedLengths(clusterZ{k}, predictedState);
        assignedZ = mgpGen3.assignMgps(clusterZ{k}, predictedState);
        
        [mgpHandles] = mgpGen3.constructMGPs(corner, predictedState, wViewed, lViewed);

    % Update UKF
    assignedZo = reshape(assignedZ', 14,1);
    %ukf.setNrMGPS(length(assignedZ));
    ukf.updateMoments(mgpHandles, assignedZo);


    % Let's see how we did
    if k == 1
        fig = figure; fig.Position = [100 100 1200 800];
    end
    hold off
    plot(clusterZ{k}(:,1), clusterZ{k}(:,2),'rx');hold on; axis equal
    drawMyRide(predictedState, 'c')
    drawMyRide(ukf.upSt,'b')

    plot(assignedZ(:,1), assignedZ(:,2),'g*','MarkerSize',10)

    now = 1;
    while now
        keydown = waitforbuttonpress;
        if keydown == 1
            now = 0;
        end
    end


end




