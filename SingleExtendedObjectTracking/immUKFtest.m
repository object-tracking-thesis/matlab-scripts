%% Test with IMM filter 

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
f1 = @(st) [st(1)+T*st(3)*cos(st(4));...
            st(2)+T*st(3)*sin(st(4));...
            st(3);...
            st(4) + T*st(5);...
            st(5);
            st(6);
            st(7)];

f2 = @(st) [st(1)+T*st(3)*cos(st(4));...
            st(2)+T*st(3)*sin(st(4));...
            st(3);...
            st(4);...
            st(5);
            st(6);
            st(7)];      
       
%% Define motion noise for model 1 (turning)
velCov1 = 0.5^2;    % 0.05 
phiDotCov1 = 0.65^2; % 0.05
% Define width and length noise 
wCov1 = 0.03^2; % std.dev is 1 cm
lCov1 = 0.05^2; % std.dev is 2 cm
R = 0.1^2;          % 0.3
subQ1 = diag([velCov1, phiDotCov1, wCov1, lCov1]);
gamma1 = [0 0 1 0 0 0 0 ;
         0 0 0 0 1 0 0 ;
         0 0 0 0 0 1 0 ;
         0 0 0 0 0 0 1]';
 
% MOTION COVARIANCE MATRIX 
Q1 = T*gamma1*subQ1*gamma1';

%% Define motion noise for model 2 (no turning)
velCov2 = 0.5^2;    % 0.05 
phiCov2 = 0.30^2;    % Heading 
% Define width and length noise 
wCov2 = 0.03^2; % std.dev is 1 cm
lCov2 = 0.05^2; % std.dev is 2 cm

subQ2 = diag([velCov2, phiCov2, wCov2, lCov2]);
gamma2 = [0 0 1 0 0 0 0 ;
          0 0 0 1 0 0 0 ;
          0 0 0 0 0 1 0 ;
          0 0 0 0 0 0 1]';
 
% MOTION COVARIANCE MATRIX 
Q2 = T*gamma2*subQ2*gamma2';
%%
f = f1;
run('mcSimOfP0.m') 
P0 = P;  

stateStorage = cell(1,120);

%%

for k = 1:nrIter
    if k == 1
        TPM = [0.50 0.50;
               0.50 0.50];
        
        nObsSt = 2;
        nSt = 7;
        st0 = x0;
        cov0 = P0;
        rCov= R;
        
        imm = CarIMM(nObsSt, nSt, st0, cov0, f1, f2, Q1, Q2, rCov, TPM);
        
        imm.mmPredict();
        
        predictedState1 = imm.mmPredSt(:,1);
        predictedState2 = imm.mmPredSt(:,2);
    else
        
        imm.mmPredict();
        
        predictedState1 = imm.mmPredSt(:,1);
        predictedState2 = imm.mmPredSt(:,2);
    end

    % Initialize mgpGen3 and use it to create MGPs
    if k == 1
        N = 4;
        mgpGen3 = MGPgenerator3(N);
    end
    
    [mgpHandles1, assignedZ1] = mgpGen3.generate(clusterZ{k}, predictedState1);
    
    [mgpHandles2, assignedZ2] = mgpGen3.generate(clusterZ{k}, predictedState2);


    % Update IMM
    assignedZo = reshape(assignedZ1', 2*length(assignedZ1),1);

    imm.mmUpdate({mgpHandles1, mgpHandles2}, assignedZo);
    
    upSt = imm.upSt;
    
    stateStorage{k} = upSt;
    % Let's see how we did
    if k == 1
        fig = figure; fig.Position = [100 100 1200 800];
    end
    hold off
    plot(clusterZ{k}(:,1), clusterZ{k}(:,2),'rx');hold on; axis equal; grid on
    p1 = drawMyRide(predictedState1, 'c');
    p2 = drawMyRide(predictedState2, [1 0.75 0.4]);
    drawMyRide(upSt,'b')
    
    plot(assignedZ1(:,1), assignedZ1(:,2),'g*','MarkerSize',10)
    
      title(sprintf('T = %.1f/%.1f [sec]', k/10,nrIter/10),'FontSize', 24)
      xlabel('X','FontSize',24)
      ylabel('Y','FontSize',24)
      dim = [.15 .4 .3 .3];
    
      str = {sprintf('$x = %.2f$',upSt(1)),...
             sprintf('$y = %.2f$',upSt(2)),...
             sprintf('$v = %.2f$',upSt(3)),...
             sprintf('$\\phi = %.2f$',upSt(4)),...
             sprintf('$\\dot{\\phi} = %.2f$',upSt(5)),...
             sprintf('$w$ = %.2f',upSt(6)),...
             sprintf('$l$ = %.2f',upSt(7)),...
             sprintf('Turning: %.2f', imm.mixWeights(1)),...
             sprintf('No Turn: %.2f', imm.mixWeights(2))};
      
         a = annotation('textbox',dim,'String',str,'FitBoxToText','on');
            a.Interpreter = 'Latex'; a.FontSize = 14; a.BackgroundColor = [1 1 1];
            
    legend([p1(1), p2(2)], 'Turning Model', 'No Turning Model');

%     now = 1;
%     while now
%         keydown = waitforbuttonpress;
%         if keydown == 1
%             now = 0;
%         end
%     end
   pause(0.2)


end






