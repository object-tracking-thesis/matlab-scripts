clear all
load carClustersCutOff.mat
load egoPosition.mat
% Get a frame,
carClusters = carClustersCutOff;

nrIter = length(carClusters);
clusterZ = cell(1,nrIter);


for j = 1:nrIter

    Ntg = carClusters{j}(:,1:3);
    Ntg(:,1:2) = Ntg(:,1:2) - repmat([egoPosition{1}(1), egoPosition{1}(2)],length(Ntg),1);  

    clusterZ{j} = Ntg;
end

%% 
[m1, m2, uOp] = cornerPoint(clusterZ{j});
 
      c1 = uOp(1); c2 = uOp(2);
      n1 = uOp(3); n2 = uOp(4);
     
      xc = (-n1*c1 + n2*c2);
      yc = (-n2*c1 -n1*c2);
          
      angle = ones(1,4)*atan2(n2,n1) + [0 1 2 3]*pi/2;

      x0 = [-9.5, 2.2, 5, 2.9, 0, 1.8, 4.7]';
      
for j = 1%:nrIter
   plot(clusterZ{j}(:,1),clusterZ{j}(:,2),'x'); hold on; axis equal
   drawMyRide(x0)
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
phiDotCov = 0.1^2; % 0.05
% Define width and length noise 
wCov = 0.03^2; % std.dev is 1 cm
lCov = 0.05^2; % std.dev is 2 cm
R = 0.3^2;          % 0.3
subQ = diag([velCov, phiDotCov, wCov, lCov]);
gamma = [0 0 1 0 0 0 0 ;
         0 0 0 0 1 0 0 ;
         0 0 0 0 0 1 0 ;
         0 0 0 0 0 0 1]';
 
% MOTION COVARIANCE MATRIX 
Q = T*gamma*subQ*gamma';

%x0 = [-16.96, -2.28, 4.39, 4, 0, 1.74, 4.45]';

P0 = P;  

stateStorage = cell(1,120);

%% Initialize UKF and use it to predict

for k = 1:nrIter
    if k == 1
        nMGPS = 7;
        nObsSt = 2;
        nSt = 7;

        ukf = UKF(Q,R, nObsSt, nSt, x0, P0);
        tic
        ukf.predictMoments(f);
        predictedState = ukf.predSt;
    else        
        tic
        ukf.predictMoments(f);
        predictedState = ukf.predSt;    
    end

    % Initialize mgpGen3 and use it to create MGPs
    if k == 1
        N = 2;
        mgpGen3 = MGPgenerator3(N);
    end     
%         corner = mgpGen3.getCarCorner(clusterZ{k}, predictedState);
%         [wViewed, lViewed] = mgpGen3.getViewedLengths(clusterZ{k}, predictedState);
%         assignedZ = mgpGen3.assignMgps(clusterZ{k}, predictedState);
%         
%         [mgpHandles] = mgpGen3.constructMGPs(corner, predictedState, wViewed, lViewed);
    [mgpHandles, assignedZ] = mgpGen3.generate(clusterZ{k}, predictedState);

    % Update UKF
    assignedZo = reshape(assignedZ', 2*length(assignedZ),1);
    ukf.setNrMGPS(length(assignedZ));
    ukf.updateMoments(mgpHandles, assignedZo);
    toc
    stateStorage{k} = ukf.upSt;
    % Let's see how we did
    if k == 1
        fig = figure; fig.Position = [100 100 1200 800];
    end
    hold off
    plot(clusterZ{k}(:,1), clusterZ{k}(:,2),'rx');hold on; axis equal; grid on
    drawMyRide(predictedState, 'c')
    drawMyRide(ukf.upSt,'b')
    
    plot(assignedZ(:,1), assignedZ(:,2),'g*','MarkerSize',10)
    
      title(sprintf('T = %.1f/%.1f [sec]', k/10,nrIter/10),'FontSize', 24)
      xlabel('X','FontSize',24)
      ylabel('Y','FontSize',24)
      dim = [.15 .4 .3 .3];
    
      str = {sprintf('$x = %.2f$',ukf.upSt(1)),...
             sprintf('$y = %.2f$',ukf.upSt(2)),...
             sprintf('$v = %.2f$',ukf.upSt(3)),...
             sprintf('$\\phi = %.2f$',ukf.upSt(4)),...
             sprintf('$\\dot{\\phi} = %.2f$',ukf.upSt(5)),...
             sprintf('$w$ = %.2f',ukf.upSt(6)),...
             sprintf('$l$ = %.2f',ukf.upSt(7))};
      
         a = annotation('textbox',dim,'String',str,'FitBoxToText','on');
         a.Interpreter = 'Latex';
         a.FontSize = 10;
         a.BackgroundColor = [1 1 1];

%     now = 1;
%     while now
%         keydown = waitforbuttonpress;
%         if keydown == 1
%             now = 0;
%         end
%     end
    pause(0.2)


end
%%

close all

fig = figure; fig.Position = [100 100 1000 800];

for k = 1:nrIter  
    Z = clusterZ{k};
        
    mZ = mean(Z);
    height = mZ(3);
    hold off;
    plot3(clusterZ{k}(:,1), clusterZ{k}(:,2), clusterZ{k}(:,3),'k.'); axis equal; hold on; grid on
         drawMyRide3(stateStorage{k}, height, 'r')
    p1 = drawMyRide3(stateStorage{k}, height*0.99, 'r');
         drawMyRide3(stateStorage{k}, height*1.01, 'r')
         
    
    t = text(double(mZ(1)), double(mZ(2)), double(mZ(3)+1.5),sprintf('Length = %.2f m, Width = %.2f m', [stateStorage{k}(7), stateStorage{k}(6)]));
    t.FontSize = 20;
        
    axis([-60 0 -50 10 58 68]); 
    axis([mZ(1)-5 mZ(1)+5 mZ(2)-5 mZ(2)+5 mZ(3)-2 mZ(3)+2]); 
    
    
    az = 120;
    el = 25;
    view(az, el);
        
    xlabel('X','FontSize',30)
    ylabel('Y','FontSize',30)
    
    pause(0.1)
        
end














