close all
% Load carClusters
load carClusters.mat
% Get a frame,
N  = 39;
%%
close all
f = figure;
while 1
for N = 1:30%length(carClusters)
    clf
    Ntg = carClusters{N}(:,1:3);
    
    % take out ego position
    load egoPosition.mat
    Xo = egoPosition{N}(1);
    Yo = egoPosition{N}(2);
    % Compensate for ego position (i.e. put ego in origin)
    Ntg = Ntg(:,1:3) - repmat([Xo, Yo, 0],length(Ntg),1);
        
    grid on
    hold on
    plot(0,0,'sk','MarkerFace','k','MarkerSize',12)
    plot(Ntg(:,1), Ntg(:,2),'bx')
    tic
    [m1,m2,uOp, filtNtg] = cornerPoint(Ntg, 0.4, 0.4);
    toc
    c1 = uOp(1);
    c2 = uOp(2);
    n1 = uOp(3);
    n2 = uOp(4);
    xc = (-n1*c1 + n2*c2);
    yc = (-n2*c1 -n1*c2);
    angle = ones(1,4)*atan2(n2,n1) + [0 1 2 3]*pi/2;
    figure(f)
    plot(filtNtg(:,1), filtNtg(:,2),'rx')
    for i = 1:4
        plot([0 cos(angle(i))]+[xc, xc], [0 sin(angle(i))]+[yc, yc], '-g','LineWidth',2)
    end            
    
    axis([-15 5 -15 5])
    pause(0.1)
end
end