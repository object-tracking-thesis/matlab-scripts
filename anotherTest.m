% test out RectTarget.m 

% Load one of the car clusters 

load car1.mat 
load car2.mat

car1adj;
car2adj;

%% Plot car2adj to look at initial state


fig = figure; fig.Position = [100 100 1500 800];
for k = 1%:length(car1adj)
    hold off
    try
        plot(car1adj{k}(:,1), car1adj{k}(:,2),'kx'); axis equal; grid on; hold on
        [~, ~, uOp, filtNtg] = cornerPoint(car1adj{k}, 0.2, 0.5);
        plot(filtNtg(:,1), filtNtg(:,2),'cx')
        
        c1 = uOp(1); c2 = uOp(2);
        n1 = uOp(3); n2 = uOp(4);
        
        xc = (-n1*c1 + n2*c2);
        yc = (-n2*c1 -n1*c2);
        
        angle = ones(1,4)*atan2(n2,n1) + [0 1 2 3]*pi/2;
                
        Col = {'r', 'g', 'b','y'};
        for i = 1:4 % i = 2 is correct one for N = 1
            plot(1.*[0 cos(angle(i))]+[xc, xc], 1.*[0 sin(angle(i))]+[yc, yc], '-','Color',Col{i},'LineWidth',2)
            plot(0+xc, 0+yc, 'og','LineWidth',2)
        end        
        
        plot(xc, yc, 'ms','MarkerFaceColor','m');
        
    catch me
        disp(me.message)
    end
    axis([11 40 -5 3])
    title(sprintf('k = %i/%i', [k, length(car1adj)]))
    waitforkey();
    %pause(1)
end

%% 

x0 = [14.5, 0.6, 4, 6.1, 0, 1.8, 4.5]';       

rectTarget = RectTarget;
rectTarget.init(x0, []);
fig = figure; fig.Position = [100 100 1000 900];
nrIter = 38;
for k = 1:nrIter
    hold off
    rectTarget.predict()
    rectTarget.Pd;    
    
    if k == 5
        rectTarget.gating(car1adj{k});
        rectTarget.Pd;

        
        rectTarget.calcLikelihood(car1adj{k});
        pSt1 = rectTarget.predSt1; 
        pSt2 = rectTarget.predSt2;
    
        rectTarget.update()        
    end
    
    [upSt, upCov] = rectTarget.getState();
    
     
    
    D = car1adj{k};
    if isempty(car1adj{k})
        D = [12 3];
    end
    %plot(D(:,1), D(:,2),'kx') ; hold on
    drawMyRide(upSt,'b')
    try
    drawMyRide(pSt1, 'm')
    drawMyRide(pSt2, 'r')
    catch me
    end
    
    axis equal
    %axis([11 45 -5 3])
    
    title(sprintf('k = %i/%i', [k nrIter]))
    
    waitforkey
end

































