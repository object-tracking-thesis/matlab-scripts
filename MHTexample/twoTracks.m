%% Generate multiple tracks, 2D xy
rng(36) % 20,36 is good 
T = 1;
A = [1 T 0 0;
    0 1 0 0;
    0 0 1 T;
    0 0 0 1];

mu_s = [0 0 0 0];

Q_s = 0.05.*[0 0 0 0;
             0 1 0 0;
             0 0 0 0;
             0 0 0 1];

H = [1 0 0 0;
     0 0 1 0];

mu_m = [0 0];

Q_m = 0.125.*[1 0;
              0 1];

x1_init = [1 2 1 0]';
x2_init = [10 0 10 -2]';
% Generate State & Measurement Sequence for two targets

N = 12;

x1_state = ones(4,N);
x2_state = ones(4,N);
x1_m = ones(2,N);
x2_m = ones(2,N);

for j = 1:N
    if j == 1
        x1_state(:,j) = A*x1_init + mvnrnd(mu_s,Q_s)';
        x2_state(:,j) = A*x2_init + mvnrnd(mu_s,Q_s)';
        
        x1_m(:,j) = H*x1_state(:,j) + mvnrnd(mu_m,Q_m)';
        x2_m(:,j) = H*x2_state(:,j) + mvnrnd(mu_m,Q_m)';
    else
        
        x1_state(:,j) = A*x1_state(:,j-1) + mvnrnd(mu_s,Q_s)';
        x2_state(:,j) = A*x2_state(:,j-1) + mvnrnd(mu_s,Q_s)';
        
        x1_m(:,j) = H*x1_state(:,j) + mvnrnd(mu_m,Q_m)';
        x2_m(:,j) = H*x2_state(:,j) + mvnrnd(mu_m,Q_m)';
    end
    
end

% Plot the trajectories

h1 = figure('Position',[50 50 1920*0.8 1080*0.8]);
    hold on
    plot(x1_state(1,:),x1_state(3,:),'ob')
    plot(x1_init(1),x1_init(3),'ob','MarkerFaceColor','blue')
    plot(x1_m(1,:),x1_m(2,:),'xb')

    plot(x2_state(1,:),x2_state(3,:),'or')
    plot(x2_init(1),x2_init(3),'or','MarkerFaceColor','red')
    plot(x2_m(1,:),x2_m(2,:),'xr')

    xlabel('X')
    ylabel('Y')

    X1 = [x1_init x1_state];
    X2 = [x2_init x2_state];
    for k = 1:length(X1)-1
        a = 1+k;
        [xf, yf]=ds2nfu(X1(1,k:a),X1(3,k:a));
        annotation('arrow', xf,yf)
        [xf, yf]=ds2nfu(X2(1,k:a),X2(3,k:a));
        annotation('arrow', xf,yf)
    end






