close all
% Load carClusters
load carClustersCutOff.mat
load egoPosition.mat
% Get a frame,
N  = 10;
carClusters = carClustersCutOff;
%%
close all
f = figure;
f.Position = [100 100 1000 800];
angle = ones(N,4);
%while 1
for N = 1:length(carClusters)
    clf
    Ntg = carClusters{N}(:,1:3);
    
    Ntg(:,1:2) = Ntg(:,1:2) - repmat([egoPosition{1}(1), egoPosition{1}(2)],length(Ntg),1);
    
    figure(f)
    hold on; grid on
    op = 0.8;
    plot(Ntg(:,1), Ntg(:,2),'x','Color',op.*[1 1 1]);
    EGO = [egoPosition{N}(1)-egoPosition{1}(1), egoPosition{N}(2)-egoPosition{1}(2)];
    plot(EGO(1), EGO(2),'ko');
    
    text(EGO(1)-2,EGO(2)-2, 'EGO')
    
    %tic
    [m1,m2,uOp, filtNtg] = cornerPoint(Ntg, 0.5, 0.4); % 0.3 0.5
    %toc
    c1 = uOp(1); c2 = uOp(2);
    n1 = uOp(3); n2 = uOp(4);
    xc = (-n1*c1 + n2*c2);
    yc = (-n2*c1 -n1*c2);
    angle(N,:) = ones(1,4)*atan2(n2,n1) + [0 1 2 3]*pi/2;
    
    p1=plot(filtNtg(:,1), filtNtg(:,2),'kx');
    for i = 1:4 % i = 2 is correct one for N = 1             
        plot(4.*[0 cos(angle(N,i))]+[xc, xc], 4.*[0 sin(angle(N,i))]+[yc, yc], '-g','LineWidth',2)
        plot(0+xc, 0+yc, 'og','LineWidth',2)
    end
    
    % Test to draw correct angle         
    
    if N == 1
        bestChoice = 4;
        plot(4.*[0 cos(angle(N,bestChoice))]+[xc, xc], 4.*[0 sin(angle(N,bestChoice))]+[yc, yc], '-m','LineWidth',2)
        plot(0+xc, 0+yc, 'om','LineWidth',2)
    else
        diffAngle = mod(angle(N,:),2*pi) - mod(angle(N-1,bestChoice),2*pi);
        [~, bestChoice] = min(abs(diffAngle));        
        plot(4.*[0 cos(angle(N,bestChoice))]+[xc, xc], 4.*[0 sin(angle(N,bestChoice))]+[yc, yc], '-m','LineWidth',2)
        plot(0+xc, 0+yc, 'om','LineWidth',2)
        
    end
    
    title(sprintf('PRESS SPACE TO TIMESTEP\n K = %i / %i Nr of points: %i', [N, length(carClusters), length(filtNtg)]))
    
    text(double(xc-12), double(yc), sprintf('m1: %.2f\nm2: %.2f', [m1 m2]))
    leg = legend(p1, sprintf('m1: %.2f\nm2: %.2f', [m1 m2])); leg.FontSize = 20;
    
    %axis([-30 5 -50 5])
    axis equal
    now = 1;
    while now
        keydown = waitforbuttonpress;
        if keydown == 1
            now = 0;
        end
    end
    %pause(2)
end
%end


   
   
   
   
   
   