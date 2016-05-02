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
velCov = 1^2;
phiDotCov = 0.1^2;
 
subQ = [velCov, 0;
        0,      phiDotCov];
 
gamma = [0 0 1 0 0 0 0 ;
         0 0 0 0 1 0 0]';          
 
% MOTION COVARIANCE MATRIX 
Q = T*gamma*subQ*gamma';
 
% Define base state & base covariance (frame K = 7 in carClustersCutOff)
 
st0 = [-11.25; 1.9; 8; 3.20; 1.4; 1.6; 3.5];
 
%P0  = [];
 
%%  Draw
 
load carClustersCutOff.mat
load egoPosition.mat
% Get a frame,
% K  = 7;
% carClusters = carClustersCutOff;
%  
% Ntg = carClusters{K}(:,1:3);
%     
%     Ntg(:,1:2) = Ntg(:,1:2) - repmat([egoPosition{1}(1), egoPosition{1}(2)],length(Ntg),1);
%     
%     [m1,m2,uOp, filtNtg] = cornerPoint(Ntg, 0.1, 0); % 0.3 0.5
%     %toc
%     c1 = uOp(1); c2 = uOp(2);
%     n1 = uOp(3); n2 = uOp(4);
%     xc = (-n1*c1 + n2*c2);
%     yc = (-n2*c1 -n1*c2);
%     angle(K,:) = ones(1,4)*atan2(n2,n1) + [0 1 2 3]*pi/2;
%     
%     fig = figure; fig.Position = [100 100 1000 800];
%     hold on; grid on
%     p1=plot(filtNtg(:,1), filtNtg(:,2),'kx');
%     for i = 1:4 % i = 2 is correct one for N = 1             
%         plot(4.*[0 cos(angle(K,i))]+[xc, xc], 4.*[0 sin(angle(K,i))]+[yc, yc], '-g','LineWidth',2)
%         plot(0+xc, 0+yc, 'og','LineWidth',2)
%     end
    
%% Initiate a MGPgenerator, we want N = 1;
N = 0;
MGPgen = MGPgenerator(N);
 
  
%%
%  
% vehicleCase = MGPgen.decideCase(Ntg, st0);
%  
% [mgps, jacobians] = MGPgen.constructMGPs(vehicleCase, st0);
%  
% [assignedMgps, assignedJacobians, assignedZ] = MGPgen.assignMgps(mgps, jacobians, Ntg);
%  
% %%
%  figure(fig)
%      plot(assignedMgps(:,1), assignedMgps(:,2),'ro','MarkerSize',15);%,'MarkerFaceColor','r');
%      plot(assignedZ(:,1), assignedZ(:,2),'rx');%,'MarkerFaceColor','r');
% %     
%% Write Update function    
nrIter = 14;
 
stateVector = ones(7,nrIter);
covarianceVector = ones(7,7,nrIter);
 
% P0 = rand(7,7);
% P0 = P0'*P0;
a = 0; b = 1;
myMat = a + (b-a).*rand(7,7);
P0 = cov([f(st0)'; f(f(st0))'])+ myMat'*myMat./100000;

for k = 1:nrIter    
    if k == 1
        % Prediction
        stateVector(:,k) = st0;%f(st0);
        covarianceVector(:,:,k) = P0;%J(st0)*P0*J(st0)';
    else
        stateVector(:,k) = f(st0);
        covarianceVector(:,:,k) = J(st0)*P0*J(st0)';        
    end
        % Update 
        Ntg = carClusters{7+k}(:,1:3);
        Ntg(:,1:2) = Ntg(:,1:2) - repmat([egoPosition{1}(1), egoPosition{1}(2)],length(Ntg),1);
        [~,~,~, Ntg] = cornerPoint(Ntg, 0.4, 0.5); % 0.3 0.5
        Ntg = [Ntg, ones(length(Ntg),1)];
        
        vehicleCase = MGPgen.decideCase(Ntg, stateVector(:,k));  
        vehicleCase = [1 2];
        [mgps, jacobians] = MGPgen.constructMGPs(vehicleCase, stateVector(:,k));
        [assignedMgps, assignedJacobians, assignedZ] = MGPgen.assignMgps(mgps, jacobians, Ntg);
        
        % Stack them vertically
        len = (3+N*length(vehicleCase))*2;             
        vMgps = reshape(assignedMgps.',len,1);
        vZ = reshape(assignedZ.',len,1);
        
        vJacobs = reshape(permute(assignedJacobians, [1,3,2]),[],size(assignedJacobians,2));
        
        S = vJacobs*covarianceVector(:,:,k)*vJacobs' + (0.002^2)*eye(len);
        K = covarianceVector(:,:,k)*vJacobs'*inv(S);
                        
        stateVector(:,k) = stateVector(:,k) + K*(vZ - vMgps);     
end
%% 
% Let's see how we did 
carClusters = carClustersCutOff;
fig = figure; fig.Position = [100 100 1000 800];
for k = 1:nrIter
    
    
    Ntg = carClusters{7+k}(:,1:3);
    Ntg(:,1:2) = Ntg(:,1:2) - repmat([egoPosition{1}(1), egoPosition{1}(2)],length(Ntg),1);
    [~,~,~, Ntg] = cornerPoint(Ntg, 0.5, 0.4); % 0.3 0.5
    Ntg = [Ntg, ones(length(Ntg),1)];
    
    [m1,m2,uOp, filtNtg] = cornerPoint(Ntg, 0.1, 0); % 0.3 0.5
    %toc
    c1 = uOp(1); c2 = uOp(2);
    n1 = uOp(3); n2 = uOp(4);
    xc = (-n1*c1 + n2*c2);
    yc = (-n2*c1 -n1*c2);
    angle(k,:) = ones(1,4)*atan2(n2,n1) + [0 1 2 3]*pi/2;
    
    figure(fig)    
    p1=plot(filtNtg(:,1), filtNtg(:,2),'kx');
    hold on; grid on
    for i = 1:4 % i = 2 is correct one for N = 1
        plot(4.*[0 cos(angle(k,i))]+[xc, xc], 4.*[0 sin(angle(k,i))]+[yc, yc], '-g','LineWidth',2)
        plot(0+xc, 0+yc, 'og','LineWidth',2)
        
    end
    title(sprintf('k = %i',[k]));
    axis([-18 5 -5 5]); axis equal
    
    plot(stateVector(1,k), stateVector(2,k),'r*')
    
    C1 = corner1(stateVector(:,k));
    C2 = corner2(stateVector(:,k));
    C3 = corner3(stateVector(:,k));
    C4 = corner4(stateVector(:,k));

    plot([C1(1) C2(1)], [C1(2) C2(2)],'-r') 
    plot([C2(1) C3(1)], [C2(2) C3(2)],'-b') 
    plot([C3(1) C4(1)], [C3(2) C4(2)],'-r')
    plot([C4(1) C1(1)], [C4(2) C1(2)],'-b')
    plot(stateVector(1,k), stateVector(2,k),'sm')


    dirX = [stateVector(1,k) stateVector(1,k)+cos(stateVector(4,k))];
    dirY = [stateVector(2,k) stateVector(2,k)+sin(stateVector(4,k))];
    plot(dirX, dirY,'-k')
    
    
    
    hold off
    waitforbuttonpress
    
    
end



 %% Draw a car
 

corner1 = @(st) ...
    [...
    st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4));...
    st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4))
    ];

corner2 = @(st) ... 
    [...
    st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + st(6)*cos(st(4) + pi/2);...
    st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + st(6)*sin(st(4) + pi/2)
    ];

corner3 = @(st) ...
    [...
    st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + st(6)*cos(st(4) + pi/2) + st(7)*cos(st(4));...
    st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + st(6)*sin(st(4) + pi/2) + st(7)*sin(st(4))
    ];

corner4 = @(st)...
    [...
    st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + st(7)*cos(st(4));...
    st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + st(7)*sin(st(4))
    ];


% C1 = corner1(stateVector(:,1));
% C2 = corner2(stateVector(:,1));
% C3 = corner3(stateVector(:,1));
% C4 = corner4(stateVector(:,1));
%  
%  
%  hold on
% plot([C1(1) C2(1)], [C1(2) C2(2)],'-r') 
% plot([C2(1) C3(1)], [C2(2) C3(2)],'-b') 
% plot([C3(1) C4(1)], [C3(2) C4(2)],'-r')
% plot([C4(1) C1(1)], [C4(2) C1(2)],'-b')
% plot(stateVector(1,1), stateVector(2,1),'sm')
% 
% 
% dirX = [stateVector(1,1) stateVector(1,1)+cos(stateVector(4,1))];
% dirY = [stateVector(2,1) stateVector(2,1)+sin(stateVector(4,1))];
% plot(dirX, dirY,'-k')
% axis equal
%  
 
 
 
 
    

