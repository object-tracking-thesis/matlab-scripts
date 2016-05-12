% This file contains illustrative plots for the thesis 
close all
set(0,'defaulttextinterpreter','latex')

r = .5;
phi = 0:0.01:2*pi;

xc = r*cos(phi);
yc = r*sin(phi);

n = 2;
P = [n,n]; 


fig = figure; fig.Position = [250 250 800 800];
fig.Name = 'spatialExample';
    hold on
    
    plot(xc,yc,'-k')
    plot(P(1), P(2),'ko','MarkerSize',10,'MarkerFaceColor','k')
    axis([-1 2.5 -1 2.5]); axis square
    op = 0.3;
    K = 210;
    plot([xc(K) P(1)], [yc(K) P(2)],'--k','Color',op*ones(1,3))
    K = 580;
    plot([xc(K) P(1)], [yc(K) P(2)],'--k','Color',op*ones(1,3))
    plot(0,0,'k.')        
    
    for j = 600:10:(580+240)
       h = mod(j, length(phi)); 
       
       plot(xc(h)+mvnrnd(0.1,0.05)/20, yc(h)+mvnrnd(0.1,0.05)/20, 'rx')
    end
       
    
    xlabel('$X$')
    ylabel('$Y$')
        
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
    text(n, n-0.2, sprintf('Lidar\nSensor'))
    text(0.5, -0.35, sprintf('Extended\nTarget'))
    text(0.6, 0.65, sprintf('Sensor\nView'))
    

    
    %annotation('textarrow',[0.55 0.45],[0.4 0.45],'String','y = x ')

    print(fig,fig.Name,'-depsc','-tiff','-cmyk','-r0')