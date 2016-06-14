
D = [ 0 0;
     -2 1;
     -4 3;
     -3 4;
      0 6;
      2 1];

phi = 3*pi/2;
R = [cos(phi) -sin(phi);
     sin(phi)  cos(phi)];
 
D = D*R;
fig = figure; fig.Position = [100 100 900 900];
Md = mean(D);

D = D - repmat(Md, length(D), 1);

plot(D(:,1),D(:,2)); hold on; axis equal
plot(-0.25,0.1,'bo','MarkerFaceColor','b')

fl = fill(D(:,1),D(:,2),'b');
fl.FaceAlpha = 0.1;
fl.EdgeColor = 'k';


lid = [14, 1];

pl = plot(lid(1), lid(2),'ko','MarkerFaceColor','k','MarkerSize',8);
pl

axis([-6 16 -4 4])

lines = [ 1.456, 3.162;
           14.0, 1;
           -0.5, -2.8333];
      
pl = plot(lines(:,1), lines(:,2),'k:');



partyMan = [];
%
A = [ 1.5, 3.167];
B = [ 2.5, 1.167];

AB = B-A;

for k = 1:6
   Dx = A' + AB'*0.2*(k-1) + 0.1*mvnrnd([0,0],eye(2))';
   partyMan = [partyMan Dx];
   plot(Dx(1), Dx(2), 'rx')
end

%
A = [ 2.5, 1.167];
B = [ 1.5,-0.833];

AB = B-A;

for k = 1:6
   Dx = A' + AB'*0.2*(k-1) + 0.1*mvnrnd([0,0],eye(2))';
   partyMan = [partyMan Dx];
   plot(Dx(1), Dx(2), 'rx')
end

A = [ 1.5,-0.833];
B = [-0.5, -2.833];

AB = B-A;

for k = 1:10
   Dx = A' + AB'*0.1*(k-1) + 0.1*mvnrnd([0,0],eye(2))';
   partyMan = [partyMan Dx];
   plot(Dx(1), Dx(2), 'rx')
end

Mp = mean(partyMan');

plot(Mp(1), Mp(2), 'ro','MarkerFaceColor','r')

 a = text(13.5,-.25, sprintf('Lidar\nSensor'));
 a.FontSize = 14; a.LineWidth = 1;
 a.Interpreter = 'Latex';
 
 a = text(-5,-.25, sprintf('Extended\nTarget'));
 a.FontSize = 14; a.LineWidth = 1;
 a.Interpreter = 'Latex';


a = text(3,-.25, sprintf('Sensor\nView'));
a.FontSize = 14; a.LineWidth = 1;
a.Interpreter = 'Latex';

axis off


fig.PaperPositionMode = 'auto';
path = '~/Desktop/';
print(fig,strcat(path,sprintf('ettExample')),'-dpng','-r300','-opengl');






