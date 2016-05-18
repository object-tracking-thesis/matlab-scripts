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
nMGPS = 1;
nObsSt = 1;
nSt = 2;
st0 = x0;
cov0 = P0;

ukf = UKF(Q,R, nMGPS, nObsSt, nSt, st0, cov0);

disp('k = 1')
ukf.predictMoments(f);
ukf.predSt
ukf.predCov

ukf.updateMoments({h}, Z(1))
ukf.upSt
ukf.upCov

disp('k = 2')
ukf.predictMoments(f);
ukf.predSt
ukf.predCov

ukf.updateMoments({h}, Z(2))
ukf.upSt
ukf.upCov


disp('k = 3')
ukf.predictMoments(f);
ukf.predSt
ukf.predCov

ukf.updateMoments({h}, Z(3))
ukf.upSt
ukf.upCov







