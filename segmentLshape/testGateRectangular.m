% Test out gate rectangular

path = '/Users/markocotra/Google Drive/Thesis Work AF - Object Tracking in Point Cloud/Data/carClusters/carClusters_kitti_campus_01_50';
load(path);

allTg = carClusters;

tg1 = cell(1,1);
tg2 = cell(1,2);

for j = 1:length(allTg)
    
    tg1{j} = allTg{j}{1};    
    try
        tg2{j} = allTg{j}{2};        
    catch me
        tg2{j} = [];
        
    end
    
end

%% 
fig = figure; fig.Position = [100 100 1000 800];

for N = 1:50;
    %N = 1;
    hold off
    plot(tg1{N}(:,1), tg1{N}(:,2),'kx'); axis equal; hold on
    
    M = mean(tg1{N}(:,1:2));
    plot(M(1), M(2), 'm*')
    [~, ~, uOp, filtNtg] = cornerPoint(tg1{N},0.2,0.5);
    
    plot(filtNtg(:,1), filtNtg(:,2),'cx')
    
    c1 = uOp(1); c2 = uOp(2);
    n1 = uOp(3); n2 = uOp(4);
    
    xc = (-n1*c1 + n2*c2);
    yc = (-n2*c1 -n1*c2);
    
    angle = ones(1,4)*atan2(n2,n1) + [0 1 2 3]*pi/2;
    
    plot(xc, yc, 'gs')
    
%     Col = {'r', 'g', 'm', 'y'};
%     for i = 1:4 % i = 2 is correct one for N = 1
%         plot(1.*[0 cos(angle(i))]+[xc, xc], 1.*[0 sin(angle(i))]+[yc, yc], '-','Color',Col{i},'LineWidth',2)
%         plot(0+xc, 0+yc, 'og','LineWidth',2)
%     end
    
    
    

    % test gateRectangular

    %massVec1 = [cos(6.0816), sin(6.0816)];
    %massVec2 = [cos(7.6524), sin(7.6524)];
    
    % Formulate all 4 possible L-vectors, to find the ones that are along the data

    M = mean(tg1{N}(:,1:2));
    
    v1 = [1, -n1/n2]; v1 = v1./norm(v1,2);
    v2 = [1,  n2/n1]; v2 = v2./norm(v2,2);
    v3 = -1.*v1;
    v4 = -1.*v2;

    p1 = [xc yc] + v1;
    p2 = [xc yc] + v2;
    p3 = [xc yc] + v3;
    p4 = [xc yc] + v4;

    d1 = (p1(1) - M(1))^2 + (p1(2) - M(2))^2;
    d2 = (p2(1) - M(1))^2 + (p2(2) - M(2))^2;
    d3 = (p3(1) - M(1))^2 + (p3(2) - M(2))^2;
    d4 = (p4(1) - M(1))^2 + (p4(2) - M(2))^2;

    % The massVec* are vectors which are in the direction of the point cloud points, from [xc yc]

    if d1 < d3
        massVec1 = v1;
    else
        massVec1 = v3;
    end

    if d2 < d4
        massVec2 = v2;
    else
        massVec2 = v4;
    end

    plot([xc xc+massVec1(1)], [yc yc+massVec1(2)],'-g')
    plot([xc xc+massVec2(1)], [yc yc+massVec2(2)],'-g')    
    
    cluster = tg1{N};
    gateCov = 0.2^2;
    boundedCluster = gateRectangular(cluster, [xc, yc], massVec1, massVec2, gateCov);
    
    Mb = mean(boundedCluster(:,1:2));


    figure(fig)
        plot(boundedCluster(:,1), boundedCluster(:,2), 'ro')
        plot(Mb(1), Mb(2), 'ms')

    a = waitforbuttonpress;
end










%%  Historgram plots
fig = figure; fig.Position = [100 100 1000 800];
for N = 1:50;
    disp(N)
    hold off
    histogram(tg1{N}(:,3)); hold on
    [~, ~, uOp, filtNtg] = cornerPoint(tg1{N},0.2,0.5);
    histogram(filtNtg(:,3))
    title(sprintf('Nr of measurements: %.f The meanvalue: %.2f Std.dev %.2f', [length(tg1{N}(:,3)), mean(tg1{N}(:,3)), std(tg1{N}(:,3))]))
    waitforbuttonpress

end
    










