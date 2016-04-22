close all
% Load carClusters 
load carClusters.mat
% Get a frame, 
N  = 39;
Ntg = carClusters{N}(:,1:3);

% take out ego position
load egoPosition.mat
Xo = egoPosition{N}(1);
Yo = egoPosition{N}(2);
% Compensate for ego position (i.e. put ego in origin)
Ntg = Ntg(:,1:3) - repmat([Xo, Yo, 0],length(Ntg),1);
hold on
plot(Ntg(:,1), Ntg(:,2),'bx')
axis equal

% Look at sorted points, based on Z
Ntg = sortrows(Ntg,3);
% Remove tailing points
lb = 0.15;
ub = 0.4;
N1 = round(length(Ntg)*lb); % 0.15
N2 = round(length(Ntg)*ub); % 0.4

sortNtg = Ntg(N1:end-N2,1:2);
plot(sortNtg(:,1), sortNtg(:,2),'rx')

%% Find clockwise ordering for points

% Get the geometric center point of the pointCloud
m = mean(sortNtg);
% Define perpendicular line
perpLine = [1 -m(1)/m(2)]./norm([1 -m(1)/m(2)],2);

% project points onto perp line, using the formula proj_s-on-v =
% s*dot(v,s)/dot(s,s)

[r,c ] = size(sortNtg);

top = sum(sortNtg .* repmat(perpLine,r,1),2);
bot = sum(perpLine.*perpLine);

scalingFactor = top./bot;
projPoints = repmat(scalingFactor,1,c) .* repmat(perpLine,r,1);
%
% Get angle for PC1 & PC2
x = [0 1];
phi = acos(sum(perpLine.*x)/(norm(x,2)*norm(perpLine,2)));
% Rotate PC1 points by phi
R = [cos(phi) -sin(phi); sin(phi) cos(phi)];

projPointsRotated = projPoints*R';

% Add Index before sorting them

projPointsRotated = [projPointsRotated, [1:length(projPointsRotated)]'];
[~,imax] = min(projPointsRotated(:,2));
sortedProjPointsRotated = sortrows(projPointsRotated,2);

% idxV represents the mapping from sorted to unsorted points
idxV = sortedProjPointsRotated(:,3);


% Remove mean 
%sortNtg = sortNtg - repmat(mean(sortNtg), length(sortNtg), 1);

orderedNtg = sortNtg(idxV,:);

% ==== PLOT ====
    f = figure;
    figure(f)
    hold on    
    
    plot(orderedNtg(:,1), orderedNtg(:,2),'rx')
    axis equal
    
%% Calc prob of xc,yc

[H, prH] = ISEDVPE(orderedNtg, 1);
figure
p1 = plot(H(1,:),prH,'bx');
hold on
p2 = plot(H(2,:),prH,'gx');

legend([p1, p2], 'xc', 'yc')


%%
figure 
hold on
p1 = plot3(H(1,:),H(2,:),min(prH).*ones(size(H(2,:))),'x');
plot3(H(1,:),H(2,:),prH,'--.')
xl = xlabel('X');
xl.Interpreter = 'latex';
xl.FontSize = 24;

yl = ylabel('Y');
yl.Interpreter = 'latex';
yl.FontSize = 24;

zl = zlabel('$p(x_c,y_c)$');
zl.Interpreter = 'latex';
zl.FontSize = 24;

hl = legend(p1, 'Corner Points $(x_c,y_c)$');
hl.Interpreter = 'latex';
hl.FontSize = 14;
grid on

    
%% Calculate L-shape & draw it 
uOp = ISED(orderedNtg);

%[H, prH, uOp] = ISEDVPE(orderedNtg,0.001);

mR1 = zeros(1,length(uOp));
mR2 = zeros(1,length(uOp));
for h = 1%:length(uOp);

c1 = uOp(1,h);
c2 = uOp(2,h);
n1 = uOp(3,h);
n2 = uOp(4,h);

ax = 4;
% Define the two perp lines
x1 = linspace(-ax,ax,10);
y1 = -n1/n2*x1 - c1/n2;
x2 = linspace(-0.125*ax,ax,10);
y2 = n2/n1*x2 - c2/n1;

xc = (-n1*c1 + n2*c2);
yc = (-n2*c1 -n1*c2);
 f = figure;
     figure(f)
     hold on
         plot(orderedNtg(:,1), orderedNtg(:,2),'gx')
         plot(x1,y1,'-r','LineWidth',2)
         plot(x2,y2,'-r','LineWidth',2)
         % corner point
         plot(xc, yc, 'sk', 'MarkerFace', 'k')
         axis equal%([-3 2 -3 2])


% Define each vector & normalize it 

V1 = [1 -n1/n2];
V2 = [1 n2/n1];
V1 = V1/norm(V1,2);
V2 = V2/norm(V2,2);

% project points onto perp line, using the formula proj_s-on-v =
% s*dot(v,s)/dot(s,s)
V = {V1 V2};
projPoints = cell(1,2);
projC = cell(1,2);

%plot([0 V1(1)], [0 V1(2)],'k','Color', 0.5*ones(1,3))
%plot([0 V2(1)], [0 V2(2)],'k','Color', 0.5*ones(1,3))

% Project points 
[r,c ] = size(orderedNtg);

for k = 1:length(V)
    top = sum(orderedNtg .* repmat(V{k},r,1),2);    
    
    scalingFactor = top;
    projPoints{k} = repmat(scalingFactor,1,c) .* repmat(V{k},r,1);
    
    % Do the same for (xc,yc)
    top = sum([xc yc] .* V{k});
    scalingFactor = top;
    projC{k} = scalingFactor.* V{k};

end

 figure(f)
 
     plot(projPoints{1}(:,1), projPoints{1}(:,2), 'rx')
     plot(projPoints{2}(:,1), projPoints{2}(:,2), 'bx')
     plot(projC{1}(:,1), projC{1}(:,2), 'ys','MarkerFace', 'y')
     plot(projC{2}(:,1), projC{2}(:,2), 'cs','MarkerFace', 'c')
     
     
 % Direction hypothesis 
angle = ones(1,4)*atan2(n2,n1) + [0 1 2 3]*pi/2;
figure(f)
for i = 1:4
plot([0 cos(angle(i))]+[xc, xc], [0 sin(angle(i))]+[yc, yc], '-k>','LineWidth',2)
end
    

% Look at length on each projected line 

phi1 = acos(dot(V1, [0 1])); % Get rotation angle 
phi2 = acos(dot(V2, [0 1]));

R1 = [cos(phi1) -sin(phi1); sin(phi1) cos(phi1)];
R2 = [cos(phi2) -sin(phi2); sin(phi2) cos(phi2)];

projPointsR1 = projPoints{1}*R1'; % Apply rotation 
projPointsR2 = projPoints{2}*R2';

projPointsR1 = projPointsR1(abs(projPointsR1) < 1e-5 == 0); % Remove rounding errors 
projPointsR2 = projPointsR2(abs(projPointsR2) < 1e-5 == 0);

projPointsR1 = projPointsR1 - min(projPointsR1); % Set measure from zero 
projPointsR2 = projPointsR2 - min(projPointsR2);


     figure
     plot(projPointsR1, ones(size(projPointsR1)), 'rx')
     figure
     plot(projPointsR2, ones(size(projPointsR2)),  'bx')

mR1(h) = max(projPointsR1);
mR2(h) = max(projPointsR2);
end
   
%% Plots of length probs 
% figure
% plot(mR1, prH,'xb')
% hold on
% plot(mR2, prH,'xr')


%% Remove outliers (TODO)

% outLimit = 0.2;
% 
% sr1 = diff([0 ;sort(projPointsR1)]);
% 
% i1 = find(sr1 > outLimit);




















