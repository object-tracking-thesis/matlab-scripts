% An evaluation of how well extended target tracking works when using
% MGPgenerator2, for the clustered car between frame N = 31 to N = 61 
% (where initial state information is gathered from state N = 30). 

% Most of the stuff is directly copied from textExtendedTargetTracking2.m
%%
clear all;close all;clc
 
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
velCov = 0.02^2;
phiDotCov = 0.04^2;
% Define width and length noise 
wCov = 0.01^2; % std.dev is 3 cm
lCov = 0.01^2; % std.dev is 5 cm
% (z_x,z_y) measurement covariance (LiDAR)
rCov = 0.075^2; % st.dev is in m
subQ = diag([velCov, phiDotCov, wCov, lCov]);
gamma = [0 0 1 0 0 0 0 ;
         0 0 0 0 1 0 0 ;
         0 0 0 0 0 1 0 ;
         0 0 0 0 0 0 1]';
 
% MOTION COVARIANCE MATRIX 
Q = T*gamma*subQ*gamma';
 
% Define base state & base covariance (frame N = 30 in carClustersCutOff) 
st30 = [-16.96, -2.28, 4.39, 4.04, 0.25, 1.75, 4.47]';
st31 = [-17.25, -2.57, 4.41, 4.07, 0.26, 1.75, 4.47]';

P0 = cov([st30';
          st31']);
      

%% Load carClusters
load carClustersCutOff.mat
load egoPosition.mat
% Get a frame,
carClusters = carClustersCutOff;
twoFrames = cell(1,31);

for N = 31:61    
    Ntg = carClusters{N}(:,1:3);
    Ntg(:,1:2) = Ntg(:,1:2) - repmat([egoPosition{1}(1), egoPosition{1}(2)],length(Ntg),1);
    twoFrames{1,N-30} = Ntg;
end

%% Initiate a MGPgenerator2
nrMgps = 12;
mgpGen2 = MGPgenerator2(nrMgps);


%% Write Update function    

fig = figure; fig.Position = [100 100 1200 800];
    hold on; grid on; op = 0.6;   
nrIter = 31;

for k = 1:nrIter
    
    if k == 1
    % Prediction
        predictedState = f(st30); 
        predictedCov = J(st30)*P0*J(st30)' + Q;        
    else
        predictedState = f(updatedState); 
        predictedCov = J(updatedCov)*updatedCov*J(updatedCov)' + Q;                
    end
    
    try
        clf(fig)
    catch me
    end
    hold off
    p1 = drawMyRide(predictedState,'c');
    hold on; grid on;
    p2 = plot(twoFrames{k}(:,1), twoFrames{k}(:,2),'kx');
    
    % Update
    corner = mgpGen2.getCarCorner(twoFrames{k}, predictedState(:,1));
    
    [wViewed, lViewed] = mgpGen2.getViewedLengths(twoFrames{k}, predictedState);
    [orderedMgps, orderedJacobs] = mgpGen2.constructMGPs(corner, predictedState, wViewed, lViewed);
    assignedZ = mgpGen2.assignMgps(twoFrames{k}, predictedState);
    
    % Stack them vertically
    len = 2*size(assignedZ,1); % How many rows do we have
    vMgps = reshape(orderedMgps.',len,1);
    vZ = reshape(assignedZ.',len,1);
    
    vJacobs = reshape(permute(orderedJacobs, [1,3,2]),[],size(orderedJacobs,2));
    
    S = vJacobs*predictedCov*vJacobs' + rCov*eye(len);
    K = predictedCov*vJacobs'*inv(S);
    
    updatedState = predictedState + K*(vZ - vMgps);
    updatedCov = predictedCov - K*S*K';
    
    p3 = drawMyRide(updatedState,'b');
    
    p4 = plot(orderedMgps(:,1), orderedMgps(:,2),'ro');
    p5 = plot(assignedZ(:,1), assignedZ(:,2),'rx');
    
    title(sprintf('T = %.1f/%.1f [sec]', k/10,nrIter/10),'FontSize', 24)
    xlabel('X','FontSize',24)
    ylabel('Y','FontSize',24)
    
    dim = [.15 .4 .3 .3];
    
    str = {sprintf('$x = %.2f$',updatedState(1)),...
           sprintf('$y = %.2f$',updatedState(2)),...
           sprintf('$v = %.2f$',updatedState(3)),...
           sprintf('$\\phi = %.2f$',updatedState(4)),...
           sprintf('$\\dot{\\phi} = %.2f$',updatedState(5)),...
           sprintf('$w$ = %.2f',updatedState(6)),...
           sprintf('$l$ = %.2f',updatedState(7))};
    
    a = annotation('textbox',dim,'String',str,'FitBoxToText','on');
    a.Interpreter = 'Latex';
    a.FontSize = 10;
    a.BackgroundColor = [1 1 1];    
    
    axis([-16 -14 -18 0])
    axis equal;
    
    ax = gca;
    ax.XTick = -40:40;
    ax.YTick = -30:30;

    leg = legend([p1(1) p3(1) p2 p4 p5],...
                 '$\xi_{k|k-1}$',...
                 '$\xi_{k|k}$',...
                 'LiDAR point measurements',...
                 'MGPs','Measurements mapped to MGPs');
             
    leg.Location = 'northwest';
    leg.Interpreter = 'latex';
    leg.FontSize = 8;

%    print(fig,sprintf('%.3i',k),'-depsc','-tiff');
    now = 1;
    while now
        keydown = waitforbuttonpress;
        if keydown == 1
            now = 0;
        end
    end
%     
end
















