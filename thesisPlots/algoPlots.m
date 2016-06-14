meas = load('~/Desktop/cm.mat');
meas = meas.cM;


%% 
N = 30;
A = meas{N}{1};



%%
fig = figure; fig.Position = [50 50 1600 800];
for N = 1:70    
    A = carClustersCutOff{N};
    A = meas{N}{1};
    subplot(2,1,1)
        plot3(A(:,1), A(:,2), A(:,3), 'kx'); axis equal; grid on
        title(sprintf('%i',N))
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
    subplot(2,1,2)
        plot3(A(:,1), A(:,2), A(:,3), 'kx'); axis equal; grid on
        title(sprintf('%i',N))
        xlabel('X')
        ylabel('Y')
        zlabel('Z')

    waitforkey
end

%%
figure
for N = 1:50
%N = 8;
A = meas{N}{1};


mu = mean(A(:,3));
s = std(A(:,3));

ub = mu + 0.25*s;
lb = mu - 0.5*s;

Z = A(:,3);

Idx = Z<ub & Z>lb;

Zf = Z(Idx);
% histogram(Z,8); hold on
% histogram(Zf,3)
 

B = A(Idx,:);
C = A(~Idx,:);

subplot(3,1,1)
    hold off
    plot(A(:,1), A(:,2), 'kx'); hold on
    plot(B(:,1), B(:,2), 'rx'); axis equal
subplot(3,1,2)
    hold off
    plot3(A(:,1), A(:,2), A(:,3), 'kx'); hold on
    plot3(B(:,1), B(:,2), B(:,3), 'rx'); axis equal
    zoom on
subplot(3,1,3)
    hold off
    plot(C(:,1), C(:,2), 'bx'); axis equal
    
    waitforkey
end
    
    %%
    N = 1;
    figure
    for N = 7
        A = meas{N}{1};
        mu = mean(A(:,3));
        s = std(A(:,3));
        
        ub = mu + 0.25*s;
        lb = mu - 1*s;
        
        Z = A(:,3);
        
        Idx = Z<ub & Z>lb;
        
        B = A(Idx,:);
        
        %A = meas{N}{1};
        
        [m1, m2, uOp] = cornerPoint(B);
        
        c1 = uOp(1); c2 = uOp(2);
        n1 = uOp(3); n2 = uOp(4);
        
        xc = (-n1*c1 + n2*c2);
        yc = (-n2*c1 -n1*c2);
        
        angle = ones(1,4)*atan2(n2,n1) + [0 1 2 3]*pi/2;
        hold off
        fig = figure;
        plot(A(:,1), A(:,2), 'kx'); hold on; axis equal
        
        axis([15.5 20 -1.25 1.25])
        axis off
                
        fig.PaperPositionMode = 'auto';
        path = '~/Desktop/';
        print(fig,strcat(path,sprintf('pocketNew')),'-dpng','-r300','-opengl');
%         plot(B(:,1), B(:,2), 'rx');
%         
%         plot(xc, yc, 'gs');
%         for j = 1:4
%             plot([xc, xc+cos(angle(j))], [yc, yc+sin(angle(j))],'-m','LineWidth',2)
%         end
    %waitforkey
    
    end
       
%%
        
        

    N = 1;
        
    
        A = meas{N}{1};
        mu = mean(A(:,3));
        s = std(A(:,3));
        
        ub = mu + 0.25*s;
        lb = mu - 1*s;
        
        Z = A(:,3);
        
        Idx = Z<ub & Z>lb;
        
        B = A(Idx,:);
        
        %A = meas{N}{1};
        
        [m1, m2, uOp] = cornerPoint(B);
        
        c1 = uOp(1); c2 = uOp(2);
        n1 = uOp(3); n2 = uOp(4);
        
        xc = (-n1*c1 + n2*c2);
        yc = (-n2*c1 -n1*c2);
        fig = figure
        angle = ones(1,4)*atan2(n2,n1) + [0 1 2 3]*pi/2;
        hold off
        plot(A(:,1), A(:,2), 'x','Color',[.6 .6 .6]); hold on; axis equal
        plot(B(:,1), B(:,2), 'rx');
        
%         plot(xc, yc, 'gs');
%         for j = 1:4
%             plot([xc, xc+5*cos(angle(j))], [yc, yc+5*sin(angle(j))],'-g')
%         end
%         
%        w = 1.6;
        %plot([xc, xc+w*cos(angle(4))], [yc, yc+w*sin(angle(4))],'->k','LineWidth',1)
        
%        w = 1.6;
%        plot([xc, xc+w*cos(angle(3))], [yc, yc+w*sin(angle(3))],'-k','LineWidth',1)

        axis([11.5 16.5 -1.5 2.5])
        
        [x,y] = ds2nfu([xc+0.02 12.86], [yc 1.794]);
    
        a = annotation('textarrow',x,y,'String',''); 
        a.LineWidth = 1; a.Color = 'k'; a.LineWidth = 2; a.HeadWidth = 14;
        t = text(14, -0.3, '$\mathbf{v}_2$'); t.Interpreter = 'Latex'; t.FontSize = 14;
        
        
        [x,y] = ds2nfu([xc+0.02 16.44], [yc -0.57]);
    
        a = annotation('textarrow',x,y,'String',''); 
        a.LineWidth = 1; a.Color = 'k'; a.LineWidth = 2; a.HeadWidth = 14;
        t = text(12.35, 1, '$\mathbf{v}_1$'); t.Interpreter = 'Latex'; t.FontSize = 14;
        set(gca,'xticklabel',[]);
        set(gca,'yticklabel',[]);
        plot(xc, yc,'k*','MarkerSize',12)
        axis off        

        t = text(12.3, 0.1, '$P_c$'); t.Interpreter = 'Latex'; t.FontSize = 14;
        
%         fig.PaperPositionMode = 'auto';
%         path = '~/Desktop/';
%         print(fig,strcat(path,sprintf('carVec')),'-dpng','-r300','-opengl');

        
%%     
    
fig = figure; fig.Position = [50 50 1600 800];
A = meas{N}{1};


mu = mean(A(:,3));
s = std(A(:,3));

ub = mu + 0.25*s;
lb = mu - 0.5*s;
set(0,'defaulttextinterpreter','latex')

Z = A(:,3);

Idx = Z<ub & Z>lb;

Zf = Z(Idx);
% histogram(Z,8); hold on
% histogram(Zf,3)
 

B = A(Idx,:);
C = A(~Idx,:);

subplot(2,2,2)
    hold off
    plot(A(:,1), A(:,2), 'x','Color',[.6 .6 .6]); hold on
    plot(B(:,1), B(:,2), 'rx'); axis equal; grid on
            [m1, m2, uOp] = cornerPoint(B);        
        c1 = uOp(1); c2 = uOp(2);
        n1 = uOp(3); n2 = uOp(4);        
        xc = (-n1*c1 + n2*c2);
        yc = (-n2*c1 -n1*c2);        
        angle = ones(1,4)*atan2(n2,n1) + [0 1 2 3]*pi/2;        

         plot(xc, yc, 'bs');
         for j = 1:4
             plot([xc, xc+4*cos(angle(j))], [yc, yc+4*sin(angle(j))],'-b','LineWidth',2)
         end
     axis([11.5 16.5 -1 2])
     zoom on
     xl = xlabel('$x$'); xl.FontSize = 16;
     yl = ylabel('$y$'); yl.FontSize = 16;
     set(gca,'xticklabel',[]);   
     set(gca,'yticklabel',[]);  
    
     

subplot(2,2,4)
    hold off
    plot3(A(:,1), A(:,2), A(:,3), 'x','Color',[.6 .6 .6]); hold on
    plot3(B(:,1), B(:,2), B(:,3), 'rx'); axis equal; grid on;

    xl = xlabel('$x$'); xl.FontSize = 16;
    yl = ylabel('$y$'); yl.FontSize = 16;
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
    
    zl = zlabel('$z$'); zl.FontSize = 16;
    set(gca,'zticklabel',[]);
    
    
    zoom on
subplot(2,2,1)
    plot(A(:,1), A(:,2), 'x','Color',[.6 .6 .6]); hold on; axis equal; grid on
            [m1, m2, uOp] = cornerPoint(A);        
        c1 = uOp(1); c2 = uOp(2);
        n1 = uOp(3); n2 = uOp(4);        
        xc = (-n1*c1 + n2*c2);
        yc = (-n2*c1 -n1*c2);        
        angle = ones(1,4)*atan2(n2,n1) + [0 1 2 3]*pi/2;        

         plot(xc, yc, 'bs');
         for j = 1:4
             plot([xc, xc+4*cos(angle(j))], [yc, yc+4*sin(angle(j))],'-b','LineWidth',2)
         end
         axis([11.5 16.5 -1 2])
         
         xl = xlabel('$x$'); xl.FontSize = 16;
         yl = ylabel('$y$'); yl.FontSize = 16;
         set(gca,'xticklabel',[]);
         set(gca,'yticklabel',[]);


zoom on

subplot(2,2,3)
    plot3(A(:,1), A(:,2), A(:,3), 'x','Color',[.6 .6 .6]); hold on; axis equal; grid on

         xl = xlabel('$x$'); xl.FontSize = 16;
     yl = ylabel('$y$'); yl.FontSize = 16;
     set(gca,'xticklabel',[]);   
     set(gca,'yticklabel',[]);            
     
     zl = zlabel('$z$'); zl.FontSize = 16;
     set(gca,'zticklabel',[]);        
    
    zoom on



    
    

    
 
 