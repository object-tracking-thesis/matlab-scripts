%% load data
load data/kitti_campus_01_186_clusters.mat
load data/isCarMat_kitti_campus_01_186.mat
n = length(clusters);
start_seq = 1;
end_seq = 186;

%% plot the clusters
% figure
% for i=start_seq:end_seq
%     for j=1:length(clusters{i})
%         pcshow(pointCloud(clusters{i}{j}))
%         hold on
%         axis([-10 50 -50 20 -1 2])
%         zoom(2)
%     end
%     hold off
%     pause(0.5)
% end

%% simplify to 2d xy and use only those that were found to be relevant objects
clusters_xy = cell(1,n);
for i=start_seq:end_seq    
    clusters_xy{i} = [];
    for j=1:length(clusters{i})
        if isCarMat(i,j) == 3 || isCarMat(i,j) == 4 || isCarMat(i,j) == 5
            clusters_xy{i} = [clusters_xy{i} giwMeasComp(clusters{i}{j}(:,1:2),2)];
        end
    end
end

%% birth components
means = cell(1,1);
means{1} = [0 -20 2 2 0 0]';
means{2} = [-3 -16 2 2 0 0]';
means{3} = [8 10 -2 0 0 0]';
means{4} = [28 10 -2 0 0 0]';
means{5} = [22 10 2 0 0 0]';
means{6} = [18 9 -1 0 0 0]';
means{7} = [8 0 2 2 0 0]';
means{8} = [-10 -25 4 4 0 0]';
means{9} = [30 3 -6 0 0 0]';
% means{2} = [5 0 0 10 0 0]';
% means{3} = [13 0 0 0 0 0]';
% means{4} = [20 0 0 10 0 0]';
% means{5} = [18 0 0 0 0 0]';
% means{6} = [20 0 0 -15 0 0]';

%% run filter recursion
phd_filter = NNGIWETTIMMPHDfilter;
phd_filter.set_birth_rfs(means);
targets = [];
for i=start_seq:end_seq   
    i
    figure(1);
    for j=1:length(clusters_xy{i})
        color = [0.9, 0.9, 0.9];
        plot(clusters_xy{i}(j).points(1,:),clusters_xy{i}(j).points(2,:),'x','Color',color)
        axis([-10 50 -50 20])
        hold on
    end    
    phd_filter.predict;    
    phd_filter.update(clusters_xy{i});
    phd_filter.get_number_of_targets;
    est = phd_filter.get_best_estimates;
    
    %plot the target ellipse and the measured points
    for j=1:length(est)
        rng(est(j).index)
        color =  [rand, rand, rand]; 
        [mu, P, v, V] = est(j).getState();
        cov = iwishrnd(V, v);
        [x1,x2,x3] = threeSigmaOverGrid(mu(1:2),cov);                
        plot(x3(1,:),x3(2,:),' --k','Color',color)             
        testtxt = strcat('\leftarrow i: ', num2str(est(j).index), ', w: ', num2str(est(j).weight));
        text(double(mu(1)), double(mu(2)), testtxt)        
    end
    text(25, 15, strcat('frame: ', num2str(i)))
    
    hold off
    pause(0.5)
end