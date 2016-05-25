% test out MGPgenerator4 

% Load one of the car clusters 

load car1.mat 
load car2.mat

car1adj;
car2adj;


%% Plot car2adj to look at initial state


fig = figure; fig.Position = [100 100 1000 800];
for k = 1:1%:length(car2adj)
    hold off
    plot(car2adj{k}(:,1), car2adj{k}(:,2),'kx'); axis equal; grid on; hold on
    axis([15 22 -20 0])
    waitforbuttonpress
    %pause(1)
end

[m1, m2, uOp] = cornerPoint(car2adj{k}, 0.2, 0.5);

c1 = uOp(1); c2 = uOp(2);
n1 = uOp(3); n2 = uOp(4);

xc = (-n1*c1 + n2*c2);
yc = (-n2*c1 -n1*c2);

angle = ones(1,4)*atan2(n2,n1) + [0 1 2 3]*pi/2;

plot(xc,yc,'rs','MarkerFaceColor','r')

Col = {'r', 'c', 'b','y'};
for i = 1:4 % i = 2 is correct one for N = 1
    plot(1.*[0 cos(angle(i))]+[xc, xc], 1.*[0 sin(angle(i))]+[yc, yc], '-','Color',Col{i},'LineWidth',2)
    plot(0+xc, 0+yc, 'og','LineWidth',2)
end


T = 0.1; % sample time 
f = @(st) [st(1)+T*st(3)*cos(st(4));...
           st(2)+T*st(3)*sin(st(4));...
           st(3);...
           st(4) + T*st(5);...
           st(5);
           st(6);
           st(7)];

% Initial state before first frame       
x0 = [18.3, -14.5, 4.4, 1.51, 0, 1.8, 4.5]';       

drawMyRide(f(x0),'m')

run('mcSimOfP0.m')       
       
% Define motion noise
velCov = 0.5^2;    % 0.05 
phiDotCov = 0.1^2; % 0.05
% Define width and length noise 
wCov = 0.03^2; % std.dev is 1 cm
lCov = 0.05^2; % std.dev is 2 cm
R = 0.2^2;          % 0.3
subQ = diag([velCov, phiDotCov, wCov, lCov]);
gamma = [0 0 1 0 0 0 0 ;
         0 0 0 0 1 0 0 ;
         0 0 0 0 0 1 0 ;
         0 0 0 0 0 0 1]';

% MOTION COVARIANCE MATRIX 
Q = T*gamma*subQ*gamma';

%x0 = [-16.96, -2.28, 4.39, 4, 0, 1.74, 4.45]';

P0 = P;  

stateStorage = cell(1,31);

%% Initialize UKF and use it to predict
fig = figure;
for k = 1:30%:nrIter
    if k == 1
        nObsSt = 2;
        nSt = 7;

        ukf = UKF(Q,R, nObsSt, nSt, x0, P0);
        
        ukf.predictMoments(f);
        predictedState = ukf.predSt;
    else        
        
        ukf.predictMoments(f);
        predictedState = ukf.predSt;    
    end

    % Initialize mgpGen3 and use it to create MGPs
    if k == 1
        N = 1;
        Rg = 0.1^2;
        mgpGen4 = MGPgenerator4(N,Rg);
    end     
    
    lb = 0.2;
    ub = 0.5;
    gateCov = 0.2^2;
    
    clusterZ = car2adj{k};
    
             [filtClust, uOp] = mgpGen4.preFilter(clusterZ, lb, ub, gateCov);
                       corner = mgpGen4.getCarCorner(clusterZ, uOp, predictedState)  % This one needs a lot of changing 
           [wViewed, lViewed] = mgpGen4.getViewedLengths(filtClust, uOp, predictedState);
                    corner = 3
                    assignedZ = mgpGen4.assignMgpsGating(filtClust, predictedState, wViewed, lViewed);
                 [mgpHandles] = mgpGen4.constructMGPs(corner, predictedState, wViewed, lViewed);

    [gMgpHandles, gAssignedZ] = mgpGen4.checkGates(assignedZ, mgpHandles);


     
      c1 = uOp(1); c2 = uOp(2);
      n1 = uOp(3); n2 = uOp(4);
     
      xc = (-n1*c1 + n2*c2);
      yc = (-n2*c1 -n1*c2);
    
    % Update UKF
     ukf.upSt(6)
     ukf.upSt(7)
     assignedZo = reshape(gAssignedZ', 2*length(gAssignedZ),1);          
     
     ukf.updateMoments(gMgpHandles, assignedZo);
     
     % Let's see how well we did
     figure(fig);
     hold off     
        plot(clusterZ(:,1), clusterZ(:,2), 'kx')   
     hold on; axis equal

        p1 = plot(gAssignedZ(:,1), gAssignedZ(:,2),'rx'); 

        p2 = plot(ukf.yPred(1,:), ukf.yPred(2,:),'ms');        
        
        plot(xc, yc, 'gs', 'MarkerFaceColor','g')

     drawMyRide(predictedState, 'c')
     
     drawMyRide(ukf.upSt)
%     
%     stateStorage{k} = ukf.upSt;
     waitforkey

end











