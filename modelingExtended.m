% Look at modelling the density of an extended target
figure('Position',[600 250 1000 800])
L = 4.58; % meters
B = 1.80; % meters
%


means = zeros(2,4);
covs  = zeros(2,2,4);

% Place corner gaussians

sbl = [0 0]'; % L bottom left corner
Pbl = 0.01*[0.1 0;
    -0.1 0.1];

sbr = [B 0]'; % L bottom right corner
Pbr = 0.01*[0.1 0.1;
    0   0.1];

sul = [0 L]'; % L upper left corner
Pul = 0.01*[0.1 0.1;
    0   0.1];

sur = [B L]'; % L upper right corner
Pur = 0.01*[0.1 0;
    -0.1 0.1];

means(:,1) = sbl;
means(:,2) = sbr;
means(:,3) = sul;
means(:,4) = sur;

covs(:,:,1) = Pbl;
covs(:,:,2) = Pbr;
covs(:,:,3) = Pul;
covs(:,:,4) = Pur;

hold on
n = 100; % Number of grid points
phi = linspace(0,2*pi,100); % Create grid in interval [0,2*pi] 4
% Sigma ellipse
Col = ['r', 'm','g'];

for j = 1:length(means)
    for sig = 1:3
        x = repmat(means(:,j),1,n)+sig*sqrtm(covs(:,:,j))*[cos(phi);sin(phi)]; % Equation(12) in HA1 document
        plot(x(1,:),x(2,:),Col(sig),'LineWidth',1)
    end
end
K = 30;
lDiff = L/(K+1);
bDiff = B/(K/2+1);

Lmeans = zeros(2,K);
Lcovs = zeros(2,2,K);
Bmeans = zeros(2,K);
Bcovs = zeros(2,2,K);

for k = 1:K/2
    sb = [bDiff*k 0]'; % L bottom left corner    
    Pb = 0.01*[0.1 0;
               0   0.1];
    Bmeans(:,k) = sb;
    Bcovs(:,:,k) = Pb;
    
           
    hold on;grid on
    n = 100; % Number of grid points
    phi = linspace(0,2*pi,100); % Create grid in interval [0,2*pi] 4
    % Sigma ellipse
    Col = ['r', 'm','g'];
    
    for j = 1:length(means)
        for sig = 1:3
            x = repmat(sb,1,n)+sig*sqrtm(Pb)*[cos(phi);sin(phi)]; % Equation(12) in HA1 document
            plot(x(1,:),x(2,:),Col(sig),'LineWidth',1)            
        end
    end
    
end

for k = 1:K
    sl = [0 lDiff*k]'; % L bottom left corner    
    Pl = 0.01*[0.1 0;
               0   0.1];
    
    Lmeans(:,k) = sl;
    Lcovs(:,:,k) = Pl; 
    
    hold on
    n = 100; % Number of grid points
    phi = linspace(0,2*pi,100); % Create grid in interval [0,2*pi] 4
    % Sigma ellipse
    Col = ['r', 'm','g'];
    
    for j = 1:length(means)
        for sig = 1:3
            x = repmat(sl,1,n)+sig*sqrtm(Pl)*[cos(phi);sin(phi)]; % Equation(12) in HA1 document
            plot(x(1,:),x(2,:),Col(sig),'LineWidth',1)
        end
    end
end

X = [0 0 B B 0;
    0 L L 0 0]';
plot(X(:,1),X(:,2),'b','LineWidth',2)
axis equal
%axis([-2 8 -2 8])


%% Try to rotate 

phi = angle(4)%pi/4;

R = [cos(phi) -sin(phi); sin(phi) cos(phi)];


for k = 1:K
    
    Lmeans(:,k) = R*Lmeans(:,k) + [xc;yc];
    %Lcovs(:,:,k) = R*Lcovs(:,:,k);
    
    Bmeans(:,k) = R*Bmeans(:,k) + [xc;yc];
    %Bcovs(:,:,k) = R*Bcovs(:,:,k);
end


%%

figure(f)
hold on

for k = 1:K

    n = 100; % Number of grid points
    phi = linspace(0,2*pi,100); % Create grid in interval [0,2*pi] 4
    % Sigma ellipse
    Col = ['r', 'm','c'];
    
    for j = 1:length(means)
        for sig = 1:3
            x = repmat(Lmeans(:,k),1,n)+sig*sqrtm(Lcovs(:,:,k))*[cos(phi);sin(phi)]; % Equation(12) in HA1 document
            plot(x(1,:),x(2,:),Col(sig),'LineWidth',1)
        end
    end
end


for k = 1:K/2

    hold on;grid on
    n = 100; % Number of grid points
    phi = linspace(0,2*pi,100); % Create grid in interval [0,2*pi] 4
    % Sigma ellipse
    Col = ['r', 'm','c'];
    
    for j = 1:length(means)
        for sig = 1:3
            x = repmat(Bmeans(:,k),1,n)+sig*sqrtm(Bcovs(:,:,k))*[cos(phi);sin(phi)]; % Equation(12) in HA1 document
            plot(x(1,:),x(2,:),Col(sig),'LineWidth',1)            
        end
    end
    
end




%%

