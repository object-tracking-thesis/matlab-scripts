% Load Raw Lidar Data 

pth = '~/Google Drive/Thesis Work AF - Object Tracking in Point Cloud/Data/kitti_campus_frames_rotated.mat';

data = load(pth);
%
data = data.lidarFrames;
%%

pth = '/Users/markocotra/thesis/matlab-scripts/data/kitti_campus_ellip_phd_estimates.mat';

ellip = load(pth);

ellip = ellip.ellip_phd_estimates;

tg = ellip{1}(1);

scale = tg.scale;
dof = tg.dof;

P = iwishrnd(scale, dof);

[x1, x2, x3] = threeSigmaOverGrid(tg.state(1:2,1), P);

plot(x1(1,:), x1(2,:)); hold on; axis equal
plot(x2(1,:), x2(2,:))
plot(x3(1,:), x3(2,:))



%%
confC = load('~/Desktop/confirmedTargets.mat');
confC = confC.confirmedTargets;

close all

fig = figure; fig.Position = [50 50 2560 1600];
N = length(clusters);
k = 140;
for k = 1:N
   hold off
   
   % Draw all measurements 
   plot3(data{k}(:,1), data{k}(:,2), data{k}(:,3),'.','Color',[0.7 0.7 0.7],'MarkerSize',5); hold on; axis equal
   
   % Draw car cubes 
   for h = 1:length(confC{k})
       drawMyRideCube(confC{k}{h}, -1.7, 0.7,'r')
   end
         
   % Draw ellips   
   for i = 1:length(ellip{k})
       tg = ellip{k}(i);       
       P = iwishrnd(tg.scale, tg.dof);       
       [x1, x2, x3] = threeSigmaOverGrid(tg.state(1:2,1), P);       
       Z = -1.7*ones(size(x3()));
       plot3(x3(1,:), x3(2,:), Z(1,:));
       f = fill3(x3(1,:), x3(2,:), Z(1,:), '-r');
           f.EdgeColor = 'r';
           f.FaceAlpha = 0.1;
   end
      
   % Draw all clusters
   for h = 1:length(clusters{k}) 
       p = plot3(clusters{k}{h}(:,1), clusters{k}{h}(:,2), clusters{k}{h}(:,3),'x','MarkerSize',1);
       
       % Set color depending on type 
       type = isCarMat(k,h);
       typeColors = {'k', 'g', [1 .5 0], 'c', 'b'};
       p.Color = typeColors{type};
       
   end
   
   % Draw ellips state
   for i = 1:length(ellip{k})
       tg = ellip{k}(i);              
       st = tg.state;
       p = plot3(double(st(1))*[1 1], double(st(2))*[1 1], [-1.7 2.5], 'r');
       p.LineWidth = 1;       
        txt = sprintf('id: %i', [tg.index]);
        t = text(double(st(1)), double(st(2)), 2.7, txt);              
        t.FontSize = 14;
        t.HorizontalAlignment = 'center';
        t.Interpreter = 'Latex';
        t.EdgeColor = 'r'; t.Color = 'r';
   end
%    
      % Draw car state
   for h = 1:length(confC{k})       
       st = confC{k}{h};
       
       p = plot3(double(st(1))*[1 1], double(st(2))*[1 1], [-1.7 1.9], 'r');
       p.LineWidth = 1;       
        txt = sprintf('    v: %2.f km/h \n w: %.2fm l: %.2fm', [3.6*double(st(3)) double(st(6)) double(st(7))]);
        txt = sprintf('%2.f km/h\n%.2fx%.2f m', [3.6*double(st(3)) double(st(6)) double(st(7))]);
        t = text(double(st(1)), double(st(2)), 3, txt);              
        t.FontSize = 14;
        t.HorizontalAlignment = 'center';
        t.Interpreter = 'Latex';
        t.EdgeColor = 'r'; t.Color = 'r';
   end
   
   
   %% Legend plot
   p1 = plot(0,0,'.k','MarkerSize',30); p1.Visible = 'off';
   p2 = plot(0,0,'.g','MarkerSize',30); p2.Visible = 'off';
   p3 = plot(0,0,'.','Color',[1 .5 0],'MarkerSize',30); p3.Visible = 'off';
   p4 = plot(0,0,'.c','MarkerSize',30); p4.Visible = 'off';
   p5 = plot(0,0,'.b','MarkerSize',30); p5.Visible = 'off';
   p6 = plot(0,0,'.','Color',[0.7 0.7 0.7],'MarkerSize',30); p6.Visible = 'off';
   l = legend([p6 p1 p2 p3 p4 p5],'Raw Lidar Data','Clutter', 'Car', 'Cyclist', 'Pedestrian', 'Pedestrian Group');   
   l.FontSize = 20; l.Interpreter = 'Latex';
   
   %% Axis commands  
   axis([-20 50 -50 20 -2 2])
   axis([-20 40 -35 20 -2 2])
   tx = xlabel('x'); tx.FontSize = 20;
   ty = ylabel('y'); ty.FontSize = 20;
   t = title(sprintf('T: %.1f/%.1f [sec]', [k/10  N/10])); t.FontSize = 20;
   t.Position = [30 20]; %t.Interpreter = 'Latex';
   zoom(2.5)
   %grid on
   axis off
%   waitforkey
   pause(1)   
   fig.PaperPositionMode = 'auto';
   path = '~/Desktop/demoMovie/';
   print(fig,strcat(path,sprintf('%.4i',k+20)),'-dpng','-r300','-opengl');

   k
   
   
end











