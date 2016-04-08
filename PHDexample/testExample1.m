% CV model example 
% Two targets, no clutter
rng(1337)
x1_0 = [0 1 0 1]';
x2_0 = [0 1.3 0.5 0.7]';

T = 0.1;
A = [1 T 0 0;...
     0 1 0 0;...
     0 0 1 T;...
     0 0 0 1];
 
Q = 0.1 * diag([0 1 0 1]);
%%

H = [1 0 0 0;...
     0 0 1 0];
 
R = 0.00015*diag([1 1]);

N = 12;

X1 = ones(4,N+1);
X2 = ones(4,N+1);

X1(:,1) = x1_0;
X2(:,1) = x2_0;

m1 = ones(2,N);
m2 = ones(2,N);

for k = 1:N
    if k == 1
        X1(:,k+1) = A*x1_0 + mvnrnd([0 0 0 0]',Q)';
        X2(:,k+1) = A*x2_0 + mvnrnd([0 0 0 0]',Q)';
        
        m1(:,k) = H*X1(:,k+1) + mvnrnd([0 0]',R)';
        m2(:,k) = H*X2(:,k+1) + mvnrnd([0 0]',R)';
    else
        X1(:,k+1) = A*X1(:,k) + mvnrnd([0 0 0 0]',Q)';
        X2(:,k+1) = A*X2(:,k) + mvnrnd([0 0 0 0]',Q)';
        
        m1(:,k) = H*X1(:,k+1) + mvnrnd([0 0]',R)';
        m2(:,k) = H*X2(:,k+1) + mvnrnd([0 0]',R)';
    end
end
 
%
% Define birth pdf 
mu = [0 0 0.35 0]';
pp = 0.1;
P = pp*[1 0.1 0 0;...
          0.1 1 0 0;...
          0 0 1 0.1;...
          0 0 0.1 1];
%
%

h = figure;
    hold on
    plot(X1(1,1), X1(3,1),'ob','MarkerFaceColor','b')
    plot(X1(1,1:2), X1(3,1:2),'--ob')
    plot(X1(1,2:end), X1(3,2:end),'-ob')    
    plot(m1(1,:), m1(2,:),'xb')
    
    plot(X2(1,1), X2(3,1),'ob','MarkerFaceColor','r')
    plot(X2(1,1:2), X2(3,1:2),'--or')
    plot(X2(1,2:end), X2(3,2:end),'-or')
    plot(m2(1,:), m2(2,:),'xr')
    
    axis([-1 1.5 -0.75 1.5])
    
n = 100; % Number of grid points
phi = linspace(0,2*pi,100); % Create grid in interval [0,2*pi] 4
  
% 3 Sigma ellipse position
x1 = repmat(mu(1:2:3),1,n)+1*sqrtm(pp*diag([1 1]))*[cos(phi);sin(phi)];
x2 = repmat(mu(1:2:3),1,n)+2*sqrtm(pp*diag([1 1]))*[cos(phi);sin(phi)];
x3 = repmat(mu(1:2:3),1,n)+3*sqrtm(pp*diag([1 1]))*[cos(phi);sin(phi)];
plot(x1(1,:),x1(2,:),' --k')
plot(x2(1,:),x2(2,:),' --k')
plot(x3(1,:),x3(2,:),' --k')

% generate clutter 

clutterMeas = cell(1,N);
xa = -1;
xb = 1.5;
ya = -0.5;
yb = 1.5;

for k = 1:N
   poissValue = poissrnd(Model.K);
   Xc = xa + (xb - xa)*rand(1,poissValue);
   Yc = ya + (yb - ya)*rand(1,poissValue);
   clutterMeas{1,k} = [Xc;Yc];
end
allClutter = [clutterMeas{1:k-1}];
plot(allClutter(1,:),allClutter(2,:),'x','Color',[0.3 0.3 0.3])



%% Let's test the PHD filter 

% Gather measurements 
Z = cell(1,length(m1));
for k = 1:length(Z)
   Z{k} = [m1(:,k) m2(:,k) clutterMeas{k}]; 
end

% initialize a phd filter with birthRFS
weights = {2};
means = {mu};
covariances = {P};

PHD = PHDinstance(weights, means, covariances);
gaussComps = cell(1,length(m1));
%% 

for k = 1:N
    PHD.predict();
    PHD.update(Z{k})
    PHD.N;
    gaussComps{k} = PHD.getBestComp(0.4);
end

%% 
figure(h)
hold on
for m = 1:N
    for tg = 1:10
        try
            plot(gaussComps{m}(tg).m(1), gaussComps{m}(tg).m(3), 'sk','MarkerFaceColor','k')
        catch e
        end
    end
end





