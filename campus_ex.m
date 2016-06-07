%% load data
load data/kitti_campus_01_186_clusters.mat
load data/isCarMat_kitti_campus_01_186.mat
n = length(clusters);
start_seq = 1;
end_seq = 186;

%% plot the clusters
% figure
% for i=start_seq:end_seq
%     i
%     for j=1:length(clusters{i})
%         pcshow(pointCloud(clusters{i}{j}))
%         hold on
%         axis([-10 50 -50 20 -1 2])
%         zoom(2)
%     end
%     hold off
%     pause(0.3)
% end

%% create the measurement components
meas = cell(1,n);
carMeas = cell(1,n);
% pedMeas = cell(1,n);
% cycMeas = cell(1,n);
for i=start_seq:end_seq    
    meas{i} = [];
    carMeas{i} = [];
    for j=1:length(clusters{i})
        if isCarMat(i,j) == 2
            carMeas{i} = [carMeas{i} giwMeasComp(clusters{i}{j}(:,1:3),isCarMat(i,j))];
        elseif isCarMat(i,j) == 3 || isCarMat(i,j) == 4 || isCarMat(i,j) == 5
            meas{i} = [meas{i} giwMeasComp(clusters{i}{j}(:,1:3),isCarMat(i,j))];
        end
%         elseif isCarMat(i,j) == 3
%             cycMeas{i} = [cycMeas{i} giwMeasComp(clusters{i}{j}(:,1:3),isCarMat(i,j))];
%         elseif isCarMat(i,j) == 4
%             pedMeas{i} = [pedMeas{i} giwMeasComp(clusters{i}{j}(:,1:3),isCarMat(i,j))];
%         end
    end
end

%% birth components
ellipMeans = cell(1,9);
ellipMeans{1} = [0 -20 2 2 0 0]';
ellipMeans{2} = [-3 -16 2 2 0 0]';
ellipMeans{3} = [8 10 -2 0 0 0]';
ellipMeans{4} = [7 10 -2 0 0 0]';
ellipMeans{5} = [28 10 -2 0 0 0]';
ellipMeans{6} = [22 10 2 0 0 0]';
ellipMeans{7} = [18 9 -1 0 0 0]';
ellipMeans{8} = [8 0 2 2 0 0]';
ellipMeans{9} = [-10 -24 2 2 0 0]';
ellipMeans{10} = [27 3 -3 0 0 0]';

rectMeans = cell(1,1);
rectMeans{1} = [15 0 4 0 0 1.8 4.7]';

%% run filter recursion
ellip_phd = EllipPHDfilter;
ellip_phd.set_birth_rfs(ellipMeans);
% rect_phd = RectPHDfilter;
% rect_phd.set_birth_rfs(rectMeans);
targets = [];
ellip_phd_estimates = cell(1,end_seq);
for i=start_seq:end_seq   
    i
    estimates = [];
    figure(1);
    for j=1:length(meas{i})
        color = [0.9, 0.9, 0.9];
        plot(meas{i}(j).points(1,:),meas{i}(j).points(2,:),'x','Color',color)
        axis([-10 50 -50 20])
        hold on
    end    
    ellip_phd.predict;    
    ellip_phd.update(meas{i});
    ellip_est = ellip_phd.get_best_estimates;
    
%     rect_phd.predict;    
%     rect_phd.update(carMeas{i});
%     rect_phd.get_number_of_targets
%     rect_est = rect_phd.get_best_estimates;
    
    %plot the target ellipse and the measured points
    for j=1:length(ellip_est)
        rng(ellip_est(j).index)
        color =  [rand, rand, rand]; 
        [mu, P, v, V] = ellip_est(j).getState();
        est = giwPHDEstimate(mu, v, V, ellip_est(j).index);
        estimates = [estimates est];
        cov = iwishrnd(V, v);
        [x1,x2,x3] = threeSigmaOverGrid(mu(1:2),cov);                
        plot(x3(1,:),x3(2,:),' --k','Color',color)             
        testtxt = strcat('\leftarrow i: ', num2str(ellip_est(j).index), ', w: ', num2str(ellip_est(j).weight));
        text(double(mu(1)), double(mu(2)), testtxt)
    end
    ellip_phd_estimates{i} = estimates;
%     for j=1:length(rect_est)
%         [mu, P, v, V] = rect_est(j).getState();
%         drawMyRide(mu,'r')
%     end
    text(25, 15, strcat('frame: ', num2str(i)))
    axis equal
    
    hold off
    pause(0.01)
%     waitforkey
end