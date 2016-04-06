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

H = [1 0 0 0;...
     0 0 1 0];
 
R = 0.0001*diag([1 1]);

N = 3;

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
 
figure
    hold on
    plot(X1(1,1), X1(3,1),'ob','MarkerFaceColor','b')
    plot(X1(1,1:2), X1(3,1:2),'--ob')
    plot(X1(1,2:end), X1(3,2:end),'-ob')    
    plot(m1(1,:), m1(2,:),'xb')
    
    plot(X2(1,1), X2(3,1),'ob','MarkerFaceColor','r')
    plot(X2(1,1:2), X2(3,1:2),'--or')
    plot(X2(1,2:end), X2(3,2:end),'-or')
    plot(m2(1,:), m2(2,:),'xr')
    
    axis([-0.25 0.5 -0.25 1])













