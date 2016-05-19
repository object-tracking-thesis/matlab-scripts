%% load data
load data/kitti_crossing_114_clusters.mat
n = length(clusters);
start_seq = 25;
end_seq = 29;

%% plot the clusters
% figure
% for i=start_seq:end_seq
%     pcshow(pointCloud(clusters{i}))
%     axis([0 30 0 20 -5 5])
%     zoom(2)
%     pause(0.5)
% end

%% simplify to 2d xy and use only those that were found to be cars
clusters_xy = cell(1,n);
for i=start_seq:end_seq    
    clusters_xy{i} = [];
    for j=1:length(clusters{i})
        if isCarMat(i,j) == 2
            clusters_xy{i} = [clusters_xy{i} giwMeasComp(clusters{i}{j}(:,1:2),2)];
        end
    end
end

%% motion and measurement model
%constant velocity
T = 0.1;
theta = 1;
sigma = 2;
d = 2;
F = [1      T       (1/2)*T^2;
     0      1       T;
     0      0       exp(-T/theta)];
q = 1;
Q = sigma^2*(1-exp(-(2*T)/theta))*diag([0 0 1]);

H = [1 0 0];
R = 0.00015*diag(ones(1,d));

%% birth components
means = cell(1,4);
covariances = cell(1,4);
weights = cell(1,4);
dofs = cell(1,4);
scales = cell(1,4);
means{1} = [0 0 0 0 0 0]';
means{2} = [7 11 10 -5 0 0]';
means{3} = [-9 19.5 10 -5 0 0]';
means{4} = [-12 -0.4 -3 -10 0 0]';
covariances{1} = 0.1*diag(ones(1,3));
covariances{2} = 0.1*diag(ones(1,3));
covariances{3} = 0.1*diag(ones(1,3));
covariances{4} = 0.1*diag(ones(1,3));
weights{1} = 0.005;
weights{2} = 0.005;
weights{3} = 0.005;
weights{4} = 0.005;
dofs{1} = 7;
dofs{2} = 7;
dofs{3} = 7;
dofs{4} = 7;
scales{1} = diag([1 1]);
scales{2} = diag([1 1]);
scales{3} = diag([1 1]);
scales{4} = diag([1 1]);

%% run filter recursion
giwphd_filter = GIWPHDfilter;
giwphd_filter.set_birth_rfs(means, covariances, dofs, scales, weights);
giwphd_filter.set_model_parameters(F,Q,H,R,T);
giw_comps = [];
for i=start_seq:end_seq
    figure(1);
    for j=1:length(clusters_xy{i})
        color = [0.9, 0.9, 0.9];
        plot(clusters_xy{i}(j).points(1,:),clusters_xy{i}(j).points(2,:),'x','Color',color)
        axis([-25 25 -25 25])
        hold on
    end    
    giwphd_filter.predict;    
    giwphd_filter.update(clusters_xy{i});
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
        testtxt = strcat('\leftarrow i: ', num2str(est(j).index), ', w: ', num2str(est(j).weight));
        text(double(est(j).mu(1)), double(est(j).mu(2)), testtxt)        
    end
    text(25, 15, strcat('frame: ', num2str(i)))
    
    hold off
    pause(0.5)
end