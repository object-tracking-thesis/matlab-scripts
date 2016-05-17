%% load data
load data/bicycleClusters.mat
n = length(bicycleClusters);
start_seq = 15;
end_seq = 35;

%% plot the clusters
figure
for i=start_seq:end_seq
    pcshow(pointCloud(bicycleClusters{i}))
    axis([0 30 0 20 -5 5])
    zoom(2)
    pause(0.5)
end

%% simplify to 2d xy
bicycleClusters_xy = cell(1,n);
for i=start_seq:end_seq
    bicycleClusters_xy{i} = bicycleClusters{i}(:,1:2);
end

%% motion and measurement model
%constant velocity
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

%% birth component
mu = [mean(bicycleClusters_xy{start_seq}) 0 0]';
P = 0.1*diag(ones(1,4));
v = 7;
V = diag([1 1]);
weight = 1;
index = 1;
birth_comp = giwComp(mu, P, v, V, weight, index);

%% prediction