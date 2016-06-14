
% Car
X = [0,  0;
     11, 0;
     11, 4;
     9,  7;
     2,  7;
     0,  4;
     0,  0];
% Wheels
W1 = [1, 0;
      1, -1;
      2, -1;
      2,  0];
  
W2 = W1 + repmat([8 0],length(W1),1);  
% Window
Y = [ 0.5, 4;
      2.25, 6.5;
      8.75, 6.5;
     10.5, 4;
      0.5  4];

%
for j = 1:length(X)
    if X(j,1) > 5
        X(j,:) = X(j,:) - [2 0];
    end
    if X(j,2) > 0
        X(j,:) = X(j,:) - [0 1];
    end
end

for j = 1:length(Y)
    if Y(j,1) > 5
        Y(j,:) = Y(j,:) - [2 0];
    end
    if Y(j,2) > 0
        Y(j,:) = Y(j,:) - [0 1];
    end
end

for j = 1:length(W1)
   if W1(j,1) > 5
      W1(j,:) = W1(j,:) - [2 0];       
   end   
   %W1(j,:) = W1(j,:) - [0 1];
end

for j = 1:length(W1)
   if W2(j,1) > 5
      W2(j,:) = W2(j,:) - [2 0];       
   end   
end

% Lidar Stand
LS = [3.25 6;
      3.25 6.5;
      3.75 7;
      5.25 7;
      5.75 6.5;      
      5.75 6];

% Lidar
L = [4.25, 6+1;
     4.25, 7.5+1;
     4.75, 7.5+1;
     4.75, 6+1];

fig = figure; fig.Position = [50 50 1600 800];
    plot(X(:,1), X(:,2),'-k'); axis equal; hold on
    plot(W1(:,1), W1(:,2), '-k')
    plot(W2(:,1), W2(:,2), '-k')
    plot(Y(:,1), Y(:,2), '-k')

    plot(L(:,1), L(:,2), '-k')
    plot(LS(:,1), LS(:,2), '-k','LineWidth',2)
    
    axis([-2 28 -2 9])
    
% generate beams
figure(fig)

    baseComp = [4.7 8.4];

    phi = linspace(0,-pi/6,64);
    
    for p = phi
        beam = [10*cos(p), 10*sin(p)];
        baseComp = baseComp - [0 0.015];
        fB = [baseComp;
              baseComp + beam];
  
        plot(fB(:,1), fB(:,2),':k')
        plot(fB(1,1), fB(1,2),'ko','MarkerFaceColor','k','MarkerSize',2)
    end            
    
    t = text(14.9, 8.4, '$ch_1$'); t.Interpreter = 'Latex'; t.FontSize = 14;    
    t = text(13.4, 2.4, '$ch_{64}$'); t.Interpreter = 'Latex'; t.FontSize = 14;
    
    
    x = [0.2 0.26];
    y = [0.7 0.67];
    
    a = annotation('textarrow',x,y,'String','Lidar Stand ');a.Interpreter = 'Latex';
    a.FontSize = 14; a.LineWidth = 1;
    
    x = [0.23 0.286];
    y = [0.77 0.73];
    box off
    
    
    a = annotation('textarrow',x,y,'String','Lidar ');a.Interpreter = 'Latex';
    a.FontSize = 14; a.LineWidth = 1;
    axis off
    
%% Crazy Stuff


fig = figure; fig.Position = [0 0 1200 1200*0.7071];

phi = linspace(0,2*pi,3140);

hold on
for p = phi
    plot([0 20*cos(p)], [0 20*sin(p)],':k')
    
end
axis([-10 10 -10 10])
axis off






    
    
