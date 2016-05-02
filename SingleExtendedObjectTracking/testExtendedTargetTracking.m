% Extended Tracking with EKF 
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
velCov = 0.1^2;
phiDotCov = 0.1^2;
 
subQ = [velCov, 0;
        0,      phiDotCov];
 
gamma = [0 0 1 0 0 0 0 ;
         0 0 0 0 1 0 0]';          
 
% MOTION COVARIANCE MATRIX 
Q = T*gamma*subQ*gamma';
 
% Define base state & base covariance (frame K = 7 in carClustersCutOff)
 
st0 = [-9.6220; 1.9390; 1.0288; 3.1732; 0.1; 1.71; 2.81];
 
P0  = [];
 
%%  Draw
 
load carClustersCutOff.mat
load egoPosition.mat
% Get a frame,
K  = 7;
carClusters = carClustersCutOff;
 
Ntg = carClusters{K}(:,1:3);
    
    Ntg(:,1:2) = Ntg(:,1:2) - repmat([egoPosition{1}(1), egoPosition{1}(2)],length(Ntg),1);
    
    [m1,m2,uOp, filtNtg] = cornerPoint(Ntg, 0.1, 0); % 0.3 0.5
    %toc
    c1 = uOp(1); c2 = uOp(2);
    n1 = uOp(3); n2 = uOp(4);
    xc = (-n1*c1 + n2*c2);
    yc = (-n2*c1 -n1*c2);
    angle(K,:) = ones(1,4)*atan2(n2,n1) + [0 1 2 3]*pi/2;
    
    fig = figure; fig.Position = [100 100 1000 800];
    hold on; grid on
    p1=plot(filtNtg(:,1), filtNtg(:,2),'kx');
    for i = 1:4 % i = 2 is correct one for N = 1             
        plot(4.*[0 cos(angle(K,i))]+[xc, xc], 4.*[0 sin(angle(K,i))]+[yc, yc], '-g','LineWidth',2)
        plot(0+xc, 0+yc, 'og','LineWidth',2)
    end
    
%% Initiate a MGPgenerator, we want N = 1;
N = 0;
MGPgen = MGPgenerator(N);
 
st1 = st0;%f(st0);
 
%%
 
%vehicleCase = MGPgen.decideCase(Ntg, st1);
 
%[mgps, jacobians] = MGPgen.constructMGPs(vehicleCase, st1);
 
%[assignedMgps, assignedJacobians, assignedZ] = MGPgen.assignMgps(mgps, jacobians, Ntg);
 
%%
% figure(fig)
%     plot(assignedMgps(:,1), assignedMgps(:,2),'ro','MarkerSize',15);%,'MarkerFaceColor','r');
%     plot(assignedZ(:,1), assignedZ(:,2),'rx');%,'MarkerFaceColor','r');
%     
%% Write Update function    
nrIter = 1;
 
stateVector = ones(7,nrIter);
covarianceVector = ones(7,7,nrIter);
 
P0 = rand(7,7);
P0 = P0'*P0;
 
for k = 1:nrIter
    if k == 1
        % Prediction
        stateVector(:,k) = f(st0);
        covarianceVector(:,:,k) = J(st0)*P0*J(st0)';
    else
        stateVector(:,k) = f(st0);
        covarianceVector(:,:,k) = J(st0)*P0*J(st0)';        
    end
        % Update 
        Ntg = carClusters{k}(:,1:3);        
        vehicleCase = MGPgen.decideCase(Ntg, st1);
        [mgps, jacobians] = MGPgen.constructMGPs(vehicleCase, st1);
        [assignedMgps, assignedJacobians, assignedZ] = MGPgen.assignMgps(mgps, jacobians, Ntg);
        
        % Stack them vertically
        len = (3+N*length(vehicleCase))*2;
        vMgps = reshape(assignedMgps.',len,1);
    
        vJacobs = reshape(permute(assignedJacobians, [1,3,2]),[],size(assignedJacobians,2));
        
        S = vJacobs*covarianceVector(:,:,k)*vJacobs' + (0.1)*eye(len);
end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
    

