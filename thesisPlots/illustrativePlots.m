% This file contains illustrative plots for the thesis 
close all
set(0,'defaulttextinterpreter','latex')

r = .5;
phi = 0:0.01:2*pi;

xc = 0.6.*r*cos(phi-pi/4);
yc = 1.*r*sin(phi);

n = 1.5;
P = [n,n]; 


fig = figure; fig.Position = [250 250 1000 1000];
fig.Name = 'spatialExample';
    hold on

    plot(xc,yc,'-k')
    plot(P(1), P(2),'ko','MarkerSize',10,'MarkerFaceColor','k')
    axis([-1 2.5 -1 2.5]); axis square
    op = 0.3;
    K = 210;
    plot([xc(K) P(1)], [yc(K) P(2)],'--','Color',op*ones(1,3))
    K = 580;
    plot([xc(K) P(1)], [yc(K) P(2)],'--','Color',op*ones(1,3))
    plot(0,0,'k*')        
    
    Xsto = [];
    for j = 600:10:(580+240)
       h = mod(j, length(phi)); 
       X = [xc(h)+mvnrnd(0.1,0.05)/20, yc(h)+mvnrnd(0.1,0.05)/20];
       Xsto = [Xsto;X];
       plot(X(1), X(2), 'r.','MarkerSize',6)
    end
       
    Mid = mean(Xsto);
    plot(Mid(1), Mid(2), '*','Color',1*[1 0 0])
    
    xlabel('$X$')
    ylabel('$Y$')
        
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])    
    %set(gca,'visible','off')%axis off
    
    text(n, n-0.2, sprintf('Lidar\nSensor'),'FontSize',8)
    text(0.1, -0.35, sprintf('Extended\nTarget'),'FontSize',8)
    text(0.4, 0.5, sprintf('Sensor\nView'),'FontSize',8)
    

    
    %annotation('textarrow',[0.55 0.45],[0.4 0.45],'String','y = x ')
    path = '/Users/markocotra/thesis/matlab-scripts/thesisPlots/';
    print(fig,strcat(path,fig.Name),'-depsc','-tiff','-cmyk','-r0')
    
%% Rectangular example 


P = [0 0;
     0 0.4;
     1 0.4;
     1 0;
     0 0];

 
head = [[0 0];
        [0.25 0]];
    
M1 = mean(P(1:4,:));    

R = @(x)[cos(x), -sin(x); sin(x), cos(x)];

phi = -pi/4;

P = P*R(phi);

centerPoint = M1*R(phi);

head = head*R(phi);
head = head + repmat(centerPoint,2,1);

txtPoint = [0.04 0.04];

txtPoint1 = txtPoint*R(phi + pi);
txtPoint2 = txtPoint*R(phi + pi + pi/2);
txtPoint3 = txtPoint*R(phi);
txtPoint4 = txtPoint*R(phi + pi/2);


fig = figure; fig.Position = [250 250 1000 1000];
fig.Name = 'rectExample';
    
    plot(P(:,1), P(:,2),'-k'); hold on; box off
    plot(P(:,1), P(:,2), 'ok','MarkerFaceColor','k','MarkerSize',4)
    plot(head(:,1),head(:,2),'k-')
    
    
    text(P(1,1)+txtPoint1(1)-0.02,P(1,2)+txtPoint1(2)-0.01,'$(i)$')
    text(P(2,1)+txtPoint2(1)-0.08,P(2,2)+txtPoint2(2),'$(ii)$')
    text(P(3,1)+txtPoint3(1)-0.04,P(3,2)+txtPoint3(2)+0.01,'$(iii)$')
    text(P(4,1)+txtPoint4(1)-0.03,P(4,2)+txtPoint4(2),'$(iv)$')
    
    plot([centerPoint(:,1) centerPoint(:,1)+0.1], [centerPoint(:,2) centerPoint(:,2)],'-','Color',0.6*[1 1 1])
    
    for j = linspace(0.05, pi/4-0.05,20)
        plot([centerPoint(:,1)+0.07*cos(j)], [centerPoint(:,2)+0.07*sin(j)],'.','Color',0.6*[1 1 1],'MarkerSize',1)
    end
    
    
    plot(centerPoint(:,1), centerPoint(:,2),'k.')
    
    
    text(centerPoint(:,1)+0.075, centerPoint(:,2)+0.055,'$\phi$')
    text(centerPoint(:,1)-0.1, centerPoint(:,2)-0.06,'$(x,y)$')
    
    text(0.4, 0.35,'$l$')
    text(-0.2, 0.12,'$w$')
    
    
        axis equal; axis([-.5 1 -.25 1.25]); axis square
        xlabel('$X$')
        ylabel('$Y$')

        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])    

%     subplot(1,2,2)
%     
%         load carClustersCutOff.mat
%         load egoPosition.mat
%         % Get a frame,
%         N  = 10;
%         carClusters = carClustersCutOff;
%         [~, ~, ~, filtNtg] = cornerPoint(carClusters{30}, 0.4, 0.3);
%         plot(filtNtg(:,1),filtNtg(:,2),'r.'); box off
%         axis equal; axis square; axis([85 90 -34.5 -28.5])
% 
%         xlabel('$X$')
%         ylabel('$Y$')
% 
%         set(gca,'xtick',[])
%         set(gca,'xticklabel',[])
%         set(gca,'ytick',[])
%         set(gca,'yticklabel',[])    

%         
%         
%         
% 
 path = '/Users/markocotra/Desktop/';
     print(fig,strcat(path,fig.Name),'-depsc','-tiff','-cmyk','-r300')%,'-opengl')
% 
% 
% 
% 
% 








