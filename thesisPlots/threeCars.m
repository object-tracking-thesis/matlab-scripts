load carClusters.mat
load egoPosition.mat

%%
set(0,'defaulttextinterpreter','latex')

A = carClusters{1};

A = A - repmat(mean(A),length(A),1);

fig = figure; fig.Position = [100 100 1600 800];
    plot(A(:,1), A(:,2),'x','Color',[.6 .6 .6]);axis equal; hold on

x0 = [-2 0.5 0 0 0 1.7 4.5];
    p = drawMyRide(x0,'b');
    for j = 1:length(p)
        p(j).LineWidth = 2;
    end
    plot(0.25, -0.35,'ob','MarkerSize',12,'MarkerFaceColor','b')
    
x0 = [-2 0.5 0 0.05 0 1.7 4.5];
    p = drawMyRide(x0,'b','--');
    for j = 1:length(p)
        p(j).LineWidth = 2;
    end
    plot(0.29, -0.237,'sb','MarkerSize',12,'MarkerFaceColor','b')

x0 = [-.7 2.2 0 1.4 0 1.7 4.5];
    p = drawMyRide(x0,'r');
    for j = 1:length(p)
        p(j).LineWidth = 2;
    end
    plot(-0.2448, -0.1617,'or','MarkerSize',12,'MarkerFaceColor','r')
    
x0 = [-.7 2.2 0 1.45 0 1.7 4.5];
    p = drawMyRide(x0,'r','--');
    for j = 1:length(p)
        p(j).LineWidth = 2;
    end
    plot(-0.1273, -0.136,'sr','MarkerSize',12,'MarkerFaceColor','r')

    
x0 = [-1.8 1.2 0 2.7 0 1.7 4.5];
    p = drawMyRide(x0,[1 0.5 0]);
    for j = 1:length(p)
        p(j).LineWidth = 2;
    end
    plot(-0.1291, -0.5301,'o','Color',[1 .5 0],'MarkerSize',12,'MarkerFaceColor',[1 .5 0])
    
x0 = [-1.8 1.2 0 2.8 0 1.7 4.5];
    p = drawMyRide(x0,[1 0.5 0],'--');
    for j = 1:length(p)
        p(j).LineWidth = 2;
    end    
    plot(0.03526, -0.3546,'s','Color',[1 .5 0],'MarkerSize',12,'MarkerFaceColor',[1 .5 0])
    
    
    
    
    axis square; axis equal;
    axis([-5 1 -1 5])
    
    
    

xl = xlabel('$X$'); xl.FontSize = 20;
yl = ylabel('$Y$'); yl.FontSize = 20;

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])  
    
    
