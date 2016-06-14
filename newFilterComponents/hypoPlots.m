




fig = figure; fig.Position = [100 100 800 800];

N = 5;
A = carClustersCutOff{N};
A = A - repmat(mean(A), length(A),1);


[m1, m2, uOp, filtNtg] = cornerPoint(A, 0.3, 0.4);



c1 = uOp(1); c2 = uOp(2);
n1 = uOp(3); n2 = uOp(4);

xc = (-n1*c1 + n2*c2);
yc = (-n2*c1 -n1*c2);

angle = ones(1,4)*atan2(n2,n1) + [0 1 2 3]*pi/2;



A = filtNtg;
fl = fill(bX, bY,'b');axis equal; grid on; hold on
fl.FaceAlpha = 0.05;
fl = fill(rX, rY,'r');
fl.FaceAlpha = 0.05;
fl = fill(oX, oY,[1 .5 0]);
fl.FaceAlpha = 0.05;

plot(A(:,1), A(:,2), 'x','Color',[.6 .6 .6])
% Blue car 
st_blue = [-1.7 0 0 0 0 1.7 3.5]';
p1 = drawMyRide(st_blue,'b');

for j = 1:3
    p1(j).LineWidth = 1.5;
end
bX = p1(1).XData;
bY = p1(1).YData;

bC = [0.0500   -0.8500];
plot(bC(1), bC(2),'o','Color','b', 'MarkerFaceColor','b')

% Red car 
st_red = [-1.5 1 0 1.25*pi/2 0 1.7 3.5]';
p2 = drawMyRide(st_red,'r');
for j = 1:3
    p2(j).LineWidth = 1.5;
end
rX = p2(1).XData;
rY = p2(1).YData;

rC = [-0.0450   -0.2915];
plot(rC(1), rC(2),'o','Color','r', 'MarkerFaceColor','r')

% Orange car 
st_org = [-2.3 -0.3 0 1.05*pi 0 1.7 3.5]';
p3 = drawMyRide(st_org,[1 .5 0]);
for j = 1:3
    p3(j).LineWidth = 1.5;
end
oX = p3(1).XData;
oY = p3(1).YData;

oC = [-0.4386   -0.8658];
plot(oC(1), oC(2),'o','Color',[1 .5 0], 'MarkerFaceColor',[1 .5 0])


% selected measurements 

X = [-0.004417, -2.315 , 0.1367];
Y = [-0.7061  , -0.5744, 0.8843];
plot(X, Y, '*m','MarkerSize',8)

         xl = xlabel('$x$'); xl.FontSize = 16; xl.Interpreter = 'Latex';
         yl = ylabel('$y$'); yl.FontSize = 16; yl.Interpreter = 'Latex';
         set(gca,'xticklabel',[]);
         set(gca,'xtick',[]);         
         set(gca,'yticklabel',[]);
         set(gca,'ytick',[]);         

grid off

axis([-4.5 0.5 -1.75 3.25])

fig.PaperPositionMode = 'auto';
path = '~/Desktop/';
print(fig,strcat(path,sprintf('cornerHypo')),'-dpng','-r300','-opengl');

















