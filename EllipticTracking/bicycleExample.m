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
T = 0.1;
F = [1      T       (1/2)*T^2;
     0      1       T;
     0      0       exp(-T/theta)];
q = 1;
Q = sigma^2*(1-exp(-(2*T)/theta))*diag([0 0 1]);

H = [1 0 0];
R = 0.00015*diag(ones(1,d));

%% birth components
means = cell(1,1);
covariances = cell(1,1);
dofs = cell(1,1);
scales = cell(1,1);
weights = cell(1,1);
mean_birth = mean(bicycleClusters_xy{start_seq});
%TODO state order?
means{1} = [mean_birth(1) mean_birth(2) 4 -6 0 0]';
covariances{1} = 0.1*diag(ones(1,3));
dofs{1} = 7;
scales{1} = diag([1 1]);
weights{1} = 1;

%% run filter recursion
giwphd_filter = GIWPHDfilter;
giwphd_filter.set_birth_rfs(means, covariances, dofs, scales, weights);
giwphd_filter.set_model_parameters(F,Q,H,R,T);
giw_comps = [];
for i=start_seq:end_seq
    giwphd_filter.predict;
    meas = [giwMeasComp(bicycleClusters_xy{i})];
    giwphd_filter.update(meas);
    giwphd_filter.get_number_of_targets;
    est = giwphd_filter.get_best_estimates
%     if ~isempty(est)
%         for j=1:length(est)
%             rng(est(j).index+5)
%             color =  [rand, rand, rand];          
%             plot(est(j).mu(1), est(j).mu(2),'x','Color',color)
%             hold on
%             est(j).index
%         end
%     end
    giw_comps = [giw_comps est(1)];
end

%% plot the target with it's ellipse around
figure;
for i = 1:length(giw_comps)
    mu = giw_comps(i).mu;
    cov = iwishrnd(giw_comps(i).V, giw_comps(i).v);
    [x1,x2,x3] = threeSigmaOverGrid(mu(1:2),cov);                
    plot(x3(1,:),x3(2,:),' --k')
    axis([0 30 -10 10])
    zoom(2)
    hold on
    plot(bicycleClusters_xy{start_seq+i-1}(:,1),bicycleClusters_xy{start_seq+i-1}(:,2),'x')
    hold off
    pause(0.5)
end