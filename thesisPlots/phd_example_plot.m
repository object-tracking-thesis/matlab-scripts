%% motion model
rng(1337)
T = 0.1;
F = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];
q = 1;
Q = [(q*T^3)/3  0           (q*T^2)/2   0;
     0          (q*T^3)/3   0           (q*T^2)/2;
     (q*T^2)/2  0           q*T         0;
     0          (q*T^2)/2   0           q*T];

H = [1 0 0 0;...
     0 1 0 0];
R = 0.00015*diag([1 1]);
 
%create 2d true states, X=[x y dotx doty]'
n = 20;
X1 = zeros(4,n);
X10 = [1 1 0 0]';
X1(:,1) = X10;
X2 = zeros(4,n);
X20 = [1 1 0 0]';
X2(:,1) = X20;
X3 = zeros(4,n);
X30 = [1 1 0 0]';
X3(:,1) = X30;
P_0 = 0.1*diag(ones(1,4));

Z1 = zeros(2,n);
Z2 = zeros(2,n);
Z3 = zeros(2,n);
Z1(:,1) = H*X1(:,1) + mvnrnd([0 0],R)';
Z2(:,1) = H*X2(:,1) + mvnrnd([0 0],R)';
Z3(:,1) = H*X3(:,1) + mvnrnd([0 0],R)';

for i = 2:n
    X1(:,i) = F*X1(:,i-1) + mvnrnd([0 0 0 0],Q)';
    X2(:,i) = F*X2(:,i-1) + mvnrnd([0 0 0 0],Q)';
    X3(:,i) = F*X3(:,i-1) + mvnrnd([0 0 0 0],Q)';
    Z1(:,i) = H*X1(:,i) + mvnrnd([0 0],R)';
    Z2(:,i) = H*X2(:,i) + mvnrnd([0 0],R)';
    Z3(:,i) = H*X3(:,i) + mvnrnd([0 0],R)';
end

%create clutter
clutterMeas = cell(1,n);
xmin = min([X1(1,:) X2(1,:) X3(1,:)]);
xmax = max([X1(1,:) X2(1,:) X3(1,:)]);
ymin = min([X1(2,:) X2(2,:) X3(1,:)]);
ymax = max([X1(2,:) X2(2,:) X3(1,:)]);

allClutter = [];
for k = 1:n
   poissValue = poissrnd(2);
   Xc = xmin + (xmax - xmin)*rand(1,poissValue);
   Yc = ymin + (ymax - ymin)*rand(1,poissValue);
   clutterMeas{k} = [Xc;Yc];
   allClutter = [allClutter clutterMeas{k}];
end

% Gather measurements 
Z = cell(1,n);
for k = 1:length(Z)
   Z{k} = [Z1(:,k) Z2(:,k) Z3(:,k) clutterMeas{k}]; 
end

%plot the samples and clutter
figure
plot(X1(1,:),X1(2,:),'x-',X2(1,:),X2(2,:),'x-',X3(1,:),X3(2,:),'x-')
hold on
plot(allClutter(1,:),allClutter(2,:),'x','Color',[0.3 0.3 0.3])

%% init
plot_gaussians = [];
ps = 0.99;
pd = 0.98;
%birth rfs around the center where the tracks spawn
means = cell(1,1);
covariances = cell(1,1);
weights = cell(1,1);
means{1} = [1 1 0 0]';
covariances{1} = 0.1*diag(ones(1,4));
weights{1} = 3;

%% filter recursion
figure;
gaussians = [];
phd_filter = PHDfilter;
phd_filter.set_birth_rfs(means, covariances, weights);
phd_filter.set_model_parameters(F,Q,H,R);
for i=1:20
    phd_filter.predict;
    phd_filter.update(Z{i},zeros(1,size(Z{i},2)))
    phd_filter.get_number_of_targets
    est = phd_filter.get_best_estimates;
    if ~isempty(est)
        for j=1:length(est)
            rng(est(j).index+5)
            color =  [rand, rand, rand];          
            plot(est(j).mu(1), est(j).mu(2),'x','Color',color)
            hold on
            est(j).index
        end
    end
    gaussians = [gaussians est];
    pause(0.2)
end

%% plot the gaussians
%[X,Y] = meshgrid(min([xmin ymin])-1:0.05:max([xmax ymax])+1);
[X,Y] = meshgrid(-2:0.05:3);
pdf = zeros(size(X));
figure
for i=1:length(gaussians)
    %pdf = pdf+gaussians(i).weight.*normpdfOverGrid(gaussians(i).mu(1:2),gaussians(i).P(1:2,1:2),X,Y);    
    pdf = pdf+normpdfOverGrid(gaussians(i).mu(1:2),gaussians(i).P(1:2,1:2),X,Y);
end
pdf(find(pdf == 0)) = NaN;
surf(X,Y,pdf)
xlabel('x')
ylabel('y')
zlabel('PHD')