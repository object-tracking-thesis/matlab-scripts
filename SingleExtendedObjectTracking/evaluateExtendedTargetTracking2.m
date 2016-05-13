% evaluateExtendedTargetTracking2.m
% Run against the whole carClusters Sequence, mainly to see how one-sided
% measurement update works. 


%% Load carClusters
load carClustersCutOff.mat
load egoPosition.mat
% Get a frame,
carClusters = carClustersCutOff;
twoFrames = cell(1,31);

for N = 1:length(carClusters)
    Ntg = carClusters{N}(:,1:3);
    Ntg(:,1:2) = Ntg(:,1:2) - repmat([egoPosition{1}(1), egoPosition{1}(2)],length(Ntg),1);
    
    twoFrames{1,N} = Ntg;
end

fig = figure; fig.Position = [100 100 1000 800];

for N = 1:length(carClusters)
    hold on
    plot(twoFrames{N}(:,1), twoFrames{N}(:,2), '.')
    axis([-50 0 -45 5])
    % axis equal
    grid on
    
    ax = gca;
    ax.XTick = ax.XTick(1):2:ax.XTick(end);
    ax.YTick = ax.YTick(1):2:ax.YTick(end);
    title(sprintf('k=%i/%i',[N length(carClusters)]), 'FontSize',20);
    xlabel('X','FontSize',20)
    ylabel('Y','FontSize',20)
    
    now = 1;
    while now
        keydown = waitforbuttonpress;
        if keydown == 1
            now = 0;
        end
    end
    %pause(0.15)
    
end

%% Visual inspection of first 2 frames

N = 2;
[m1, m2, uOp] = cornerPoint(twoFrames{N});
sprintf('Length m1: %.3f', m1)
sprintf('Length m2: %.3f', m2)

c1 = uOp(1); c2 = uOp(2);
n1 = uOp(3); n2 = uOp(4);

xc = (-n1*c1 + n2*c2);
yc = (-n2*c1 -n1*c2);

angle = ones(1,4)*atan2(n2,n1) + [0 1 2 3]*pi/2;

Mid = mean(twoFrames{N});
fig = figure;
    hold on
    plot(twoFrames{N}(:,1), twoFrames{N}(:,2), '.')

    plot(xc, yc, 'sr')
    
    cols = [[1 0 0];[0 1 0];[0 0 1];[0 0 0]];
    for j = 1:4
       
       plot([xc xc+cos(angle(j))], [yc yc+sin(angle(j))], '-','Color',cols(j,:)) 
    end
    drawMyRide([-9.58, 2.19, 0, 2.97, 9, 1.6, 4.5])
    plot(Mid(1), Mid(2),'mx')
    
    modRect = [-7.43, 1.66; 
               -7.49, 1.27;
               -8.47, 1.41;
               -7.57, 0.67];
    
    
    axis([-10 -6 -1 3]); grid on
    axis equal

%% Set initial state as well as initial covariance 

st1 = [-9.45, 2.18, 1.21, 2.97, 0.8, 1.8, 4.5]';
st2 = [-9.58, 2.19, 1.22, 2.98, 0.9, 1.8, 4.5]';

st0 = st1; 
P0 = cov([st1';st2']);


%

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
velCov = 0.06^2;
phiDotCov = 0.06^2;
% Define width and length noise 
wCov = 0.03^2; % std.dev is 3 cm
lCov = 0.05^2; % std.dev is 5 cm
% (z_x,z_y) measurement covariance (LiDAR)
rCov = 0.07^2; % st.dev is in m
subQ = diag([velCov, phiDotCov, wCov, lCov]);
gamma = [0 0 1 0 0 0 0 ;
         0 0 0 0 1 0 0 ;
         0 0 0 0 0 1 0 ;
         0 0 0 0 0 0 1]';
 
% MOTION COVARIANCE MATRIX 
Q = T*gamma*subQ*gamma';


% Initiate a MGPgenerator2
nrMgps = 1;
mgpGen2 = MGPgenerator2(nrMgps);

% Write Update function    

fig = figure; fig.Position = [100 100 1200 800];
    hold on; grid on; op = 0.6;   
nrIter = length(carClusters);

for k = 1:nrIter
    
    if k == 1
    % Prediction
        predictedState = f(st0); 
        predictedCov = J(st0)*P0*J(st0)' + Q;        
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
    
    axis([-16 -14 -45 5])
    axis equal;
    
    %ax = gca;
    %ax.XTick = -40:40;
    %ax.YTick = -30:30;

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


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

