rng(1334)

T = 1;

f = @(x) [x(1)+T*x(2); x(2)]; % CV model
A = [1 T;
     0 1];
 
h = @(x) x(1);
H = [1 0];

Q = [0 0;
     0 1];

R = 1;

x0 = [0; 
      1];
  
P0 = [1 0;
      0 1];

N = 3;
X = zeros(2,N);
Z = zeros(1,N);

for k = 1:N
    if k == 1
        X(:,k) = f(x0) + mvnrnd([0 0],Q)';        
    else
        X(:,k) = f(X(:,k-1)) + mvnrnd([0 0],Q)';
    end        
    
    Z(:,k) = h(X(:,k)) + mvnrnd(0,R);
    
end
%
fig = figure;
plot(X(1,:),'-x'); hold on
plot(Z,'xr'); axis([0 5 0 10])

% KF filter 
xhat = zeros(2,N);
P = zeros(2,2,N);
Ppred = P;
xhatpred = xhat;

S = zeros(1,N);
K = zeros(2,N);
V = zeros(1,N);
for j = 1:N
   % Prediction Step
   if j == 1
       xhat(:,j) = A*x0;
       P(:,:,j) = A*P0*A' + Q;
       
       Ppred(:,:,j) = P(:,:,j);
       xhatpred(:,j) = xhat(:,j);
   else
       xhat(:,j) = A*xhat(:,j-1);
       P(:,:,j) = A*P(:,:,j-1)*A' + Q;
       
        Ppred(:,:,j) = P(:,:,j);
       xhatpred(:,j) = xhat(:,j);
   end
   % Update Step
   S(:,j) = H*P(:,:,j)*H' + R;
   K(:,j) = P(:,:,j)*H'*inv(S(:,j));
   V(:,j) = Z(:,j) - H*xhat(:,j);
   
   xhat(:,j) = xhat(:,j) + K(:,j)*V(:,j);
   P(:,:,j) = P(:,:,j) - K(:,j)*S(:,j)*K(:,j)';
end

figure(fig)
plot(xhat(1,:),'ko')

%%

ukf = UKF(Q,R);
%
estX = zeros(2,N);
preX = zeros(2,N);

    [SP, W] = ukf.generatePoints(x0, P0);
    [stPred, covPred] = ukf.predictMoments(f, SP, W);        
    [SP, W] = ukf.generatePoints(stPred, covPred);    
    [stUp, covUp] = ukf.updateMoments(h, SP, W, Z(1), stPred, covPred);
    stUp
    covUp
    
    [SP, W] = ukf.generatePoints(stUp, covUp);
    [stPred, covPred] = ukf.predictMoments(f, SP, W);
    [SP, W] = ukf.generatePoints(stPred, covPred);
    [stUp, covUp] = ukf.updateMoments(h, SP, W, Z(2), stPred, covPred);
    stUp
    covUp
    
    [SP, W] = ukf.generatePoints(stUp, covUp);
    [stPred, covPred] = ukf.predictMoments(f, SP, W);
    [SP, W] = ukf.generatePoints(stPred, covPred);
    [stUp, covUp] = ukf.updateMoments(h, SP, W, Z(3), stPred, covPred);
    stUp
    covUp
% 

    

%%
figure(fig)
plot(estX(1,:), 'ko')
plot(preX(1,:), 'mo')




