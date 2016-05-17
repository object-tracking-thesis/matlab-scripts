%% load data
load data/bicycleClusters.mat
n = length(bicycleClusters);
start_seq = 15;
end_seq = 35;

%% plot the clusters
% figure
% for i=start_seq:end_seq
%     pcshow(pointCloud(bicycleClusters{i}))
%     axis([0 30 0 20 -5 5])
%     zoom(2)
%     pause(0.5)
% end

%% simplify to 2d xy
bicycleClusters_xy = cell(1,n);
for i=start_seq:end_seq
    bicycleClusters_xy{i} = bicycleClusters{i}(:,1:2);
end

%% motion and measurement model
%constant velocity
theta = 1;
tau = 5;
sigma = 2;
d = 2;

T = 0.1;
F = [1      T       (1/2)*T^2;
     0      1       T;
     0      0       exp(-T/theta)];
q = 1;
Q = sigma^2*(1-exp(-(2*T)/theta))*diag([0 0 1]);

H = [1 0 0 0;...
     0 1 0 0];
R = 0.00015*diag(ones(1,d));

%% birth component
%TODO: state order [x dotx dotdotx y doty dotdoty]' ?
mean_birth = mean(bicycleClusters_xy{start_seq});
mu = [mean_birth(1) 4 0 mean_birth(2) -6 0]';
P = 0.1*diag(ones(1,3));
v = 7;
V = diag([1 1]);
weight = 1;
index = 1;
birth_comp = giwComp(mu, P, v, V, weight, index);

%% init
ps = 1.0;
pd = 1.0;
%TODO the birth components should not be run through the prediction
giw_components = [birth_comp];

%% prediction
for i = 1:length(giw_components)
    giw_components(i).weight = ps*giw_components(i).weight;
    giw_components(i).mu = kron(F,eye(d))*giw_components(i).mu;
    giw_components(i).P = F*giw_components(i).P*F' + Q;
    temp_v = giw_components(i).v;
    giw_components(i).v = exp(-T/tau)*giw_components(i).v;
    giw_components(i).V = ((giw_components(i).v - d-1)/(temp_v - d-1)) .* giw_components(i).V;
end

%% update components
meas = giwMeasComp(bicycleClusters_xy{start_seq+1})

%% plot the target with it's ellipse around
figure;
mu = giw_components(1).mu;
cov = iwishrnd(giw_components(1).V, giw_components(1).v);
[x1,x2,x3] = threeSigmaOverGrid(mu(1:2),cov);                
plot(x3(1,:),x3(2,:),' --k')
hold on
plot(bicycleClusters_xy{start_seq+1}(:,1),bicycleClusters_xy{start_seq+1}(:,2),'x')