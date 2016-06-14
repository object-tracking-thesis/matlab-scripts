syms T e theta

F = [1 T 0.5*T^2;
     0 1 T;
     0 0 e^(-T/theta)];
 
 
 A = [1 0;
      0 1];
  
  
%% 
rng(1331)
V = diag([1,1]);

v = {4, 10, 25};

X = cell(1,3);

for j = 1:3
    X{j} = V/(v{j}-2-1);
end
fig = figure; fig.Position = [0 0 1000 1000];
plotNr = 1;

for j = 1:3
   for k = 1:3
    subplot(3,3,plotNr)
        
        p2 = plotEllip([0 0]', V, [1 2]);hold on; axis equal;
            p2.Color = [.6 .6 .6]; p2.LineWidth = 1;
        p1 = plotEllip([0 0]', iwishrnd(V,v{k}), [1 2]); 
            p1.Color = 'r'; p1.LineStyle = '-';
        p3 = plot(0,0,'x'); p3.Color = 'k';
        p4 = plotEllip([0 0]', X{k}, [1 2]);
            p4.Color = 'k'; p4.LineStyle = '--';
        
        axis(1.25*[-1 1 -1 1])
        
        if j ~= 3                        
            set(gca,'xtick',[-1:0.5:1]);
            set(gca,'xticklabel',[]);    
        else
           xl = xlabel('x');
           xl.Interpreter = 'Latex';
           xl.FontSize = 18;
        end
        
        if k ~= 1
            set(gca,'ytick',[-1:0.5:1]);
            set(gca,'yticklabel',[]); 
        else
            yl = ylabel('y');
            yl.Interpreter = 'Latex';
            yl.FontSize = 18;
        end
        
        plotNr = plotNr +1;                

        pos = get(gca, 'Position');
        if j == 1
            pos(2) = 0.7;
        end
        pos(4) = 1/4;  %% This one here changes distance between vertical subplots
        if k == 1
            pos(1) = 0.1;
            pos(3) = 0.3;
        elseif k == 2
            pos(1) = 0.35;
            pos(3) = 0.3;
        elseif k == 3
            pos(1) = 0.6;
            pos(3) = 0.3;
        end
        %pos(3) = 0.3;
        set(gca, 'Position', pos)
        
        if j == 1
           t = title(sprintf('$\\mathbf{v} = %i$', v{k}));
           t.Interpreter = 'Latex';
           t.FontSize = 20;
        end
        grid on
   end      
   
end


fig.PaperPositionMode = 'auto';
path = '~/Desktop/';
print(fig,strcat(path,sprintf('iwish')),'-dpng','-r300','-opengl');


%%



 


subplot(3,3,1)    
        p2 = plotEllip([0 0]', V, [1 2]);hold on; axis equal;
            p2.Color = [.6 .6 .6]; p2.LineWidth = 1;
        p1 = plotEllip([0 0]', iwishrnd(V,v{1}), [1 2]); 
            p1.Color = 'b'; p1.LineStyle = '-';
        p3 = plot(0,0,'x'); p3.Color = 'k';
        p4 = plotEllip([0 0]', X{1}, [1 2]);
            p4.Color = 'r'; p4.LineStyle = '--';        
        axis(1.25*[-1 1 -1 1])        
subplot(3,3,2)    
        p2 = plotEllip([0 0]', V, [1 2]);hold on; axis equal;
            p2.Color = [.6 .6 .6]; p2.LineWidth = 1;
        p1 = plotEllip([0 0]', iwishrnd(V,v{2}), [1 2]); 
            p1.Color = 'b'; p1.LineStyle = '-';
        p3 = plot(0,0,'x'); p3.Color = 'k';
        p4 = plotEllip([0 0]', X{2}, [1 2]);
            p4.Color = 'r'; p4.LineStyle = '--';        
        axis(1.25*[-1 1 -1 1])
subplot(3,3,3)    
        p2 = plotEllip([0 0]', V, [1 2]);hold on; axis equal;
            p2.Color = [.6 .6 .6]; p2.LineWidth = 1;
        p1 = plotEllip([0 0]', iwishrnd(V,v{3}), [1 2]); 
            p1.Color = 'b'; p1.LineStyle = '-';
        p3 = plot(0,0,'x'); p3.Color = 'k';
        p4 = plotEllip([0 0]', X{3}, [1 2]);
            p4.Color = 'r'; p4.LineStyle = '--';
        
        axis(1.25*[-1 1 -1 1])
       
       
       



















