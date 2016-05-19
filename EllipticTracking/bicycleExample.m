%% load data
load data/bicycleClusters.mat
n = length(bicycleClusters);
start_seq = 1;
end_seq = 50;

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
    if ~isempty(bicycleClusters{i})
        bicycleClusters_xy{i} = bicycleClusters{i}(:,1:2);
    else
        bicycleClusters_xy{i} = [];
    end
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
weights{1} = 0.9;

%% run filter recursion
giwphd_filter = GIWPHDfilter;
giwphd_filter.set_birth_rfs(means, covariances, dofs, scales, weights);
giwphd_filter.set_model_parameters(F,Q,H,R,T);
giw_comps = [];
for i=start_seq:end_seq
    figure(1);
    giwphd_filter.predict;
    meas = [giwMeasComp(bicycleClusters_xy{i})];
    giwphd_filter.update(meas);
    giwphd_filter.get_number_of_targets;
    est = giwphd_filter.get_best_estimates;
    
    %plot the target ellipse and the measured points
    for j=1:length(est)
        rng(est(j).index)
        color =  [rand, rand, rand]; 
        mu = est(j).mu;
        cov = iwishrnd(est(j).V, est(j).v);
        [x1,x2,x3] = threeSigmaOverGrid(mu(1:2),cov);                
        plot(x3(1,:),x3(2,:),' --k','Color',color)
        hold on
        axis([0 30 -10 20])
        %plot(est(j).mu(1), est(j).mu(2),'x','Color',color)
        testtxt = strcat('\leftarrow i: ', num2str(est(j).index), ', w: ', num2str(est(j).weight));
        text(double(est(j).mu(1)), double(est(j).mu(2)), testtxt)
    end
    if ~isempty(bicycleClusters_xy{start_seq+i-1})
        plot(bicycleClusters_xy{start_seq+i-1}(:,1),bicycleClusters_xy{start_seq+i-1}(:,2),'x')
    end
    
    hold off
    pause(0.5)
    giw_comps = [giw_comps est];
end

%% plot the target with it's ellipse around
% figure;
% for i = 1:length(giw_comps)
%     mu = giw_comps(i).mu;
%     cov = iwishrnd(giw_comps(i).V, giw_comps(i).v);
%     [x1,x2,x3] = threeSigmaOverGrid(mu(1:2),cov);                
%     plot(x3(1,:),x3(2,:),' --k')
%     axis([0 30 -10 20])
%     zoom(1)
%     hold on
%     if ~isempty(bicycleClusters_xy{start_seq+i-1})
%         plot(bicycleClusters_xy{start_seq+i-1}(:,1),bicycleClusters_xy{start_seq+i-1}(:,2),'x')
%     end
%     hold off
%     pause(0.5)
% end