% Test to write occlusion model 
rng(1337)
% Generate random 2D point clouds 
N = 15;
PC1 = mvnrnd([4 3]', 0.1.*[1 0.2;0.2 1],N);

PC2 = mvnrnd([-2 2]', 0.1.*[1 0.2;0.2 1],N);

figure('Position',[600 100 1000 800])
    hold on
    p1 = plot(0,0,'ko','MarkerFaceColor','k');
    legend(p1,'EGO Vehicle')
    % PC1
    plot(PC1(:,1), PC1(:,2),'.','MarkerSize',20)
    m1 = mean(PC1);
    plot(m1(1),m1(2),'bx','MarkerSize',12)
    %plot([0 m1(1)], [0 m1(2)],'b')
    % define perpendicular line
    A1 = [1 -m1(1)/m1(2)]./norm([1 -m1(1)/m1(2)],2);    
    
    %plot([m1(1) m1(1)+A1(1)], [m1(2) m1(2)+A1(2)],'b')    
   
    % PC2
    plot(PC2(:,1), PC2(:,2),'.','MarkerSize',20)        
    m2 = mean(PC2);
    plot(m2(1),m2(2),'rx','MarkerSize',12)  
    %plot([0 m2(1)], [0 m2(2)],'r')
    % define perpendicular line         
    A2 = [1 -m2(1)/m2(2)]./norm([1 -m2(1)/m2(2)],2);
        
    %plot([m2(1) m2(1)+A2(1)], [m2(2) m2(2)+A2(2)],'r')      

    axis(0.5*[-10 10 -10 10])
    
    
% project points from PC2 & PC1 onto perp line
proj1 = zeros(size(PC1));
proj2 = zeros(size(PC2));

for k = 1:size(proj1,1)
    proj1(k,:) = sum(PC1(k,:).*A1) / sum(A1.*A1) .* A1;
end

for k = 1:size(proj2,1)
    proj2(k,:) = sum(PC2(k,:).*A2) / sum(A2.*A2) .* A2;
end

proj1 = repmat(m1,size(proj1,1),1) + proj1;
%plot(proj1(:,1), proj1(:,2),'db')

proj2 = repmat(m2,size(proj2,1),1) + proj2;
%plot(proj2(:,1), proj2(:,2),'dr')

% Rotate points to Y, so that we can find max & min index 

% Get angle for PC1 & PC2
x = [0 1];
phi1 = acos(sum(A1.*x)/(norm(x,2)*norm(A1,2)));
phi2 = acos(sum(A2.*x)/(norm(x,2)*norm(A2,2)));
% Rotate PC1 points by phi
R1 = [cos(phi1) -sin(phi1); sin(phi1) cos(phi1)];

R2 = [cos(phi2) -sin(phi2); sin(phi2) cos(phi2)];

proj1rot = proj1*R1';
proj2rot = proj2*R2';
% plot rotated points 
%plot(proj1rot(:,1), proj1rot(:,2),'db')
%hold on
%plot(proj2rot(:,1), proj2rot(:,2),'dr')

% Find min max points 

%[~, i1Max] = max(proj1rot(:,2))
%[~, i1Min] = min(proj1rot(:,2))


%[~, i2Max] = max(proj2rot(:,2))
%[~, i2Min] = min(proj2rot(:,2))

% plot(10.*[0 PC1(i1Max,1)],10.*[0 PC1(i1Max,2)],'m')
% plot(10.*[0 PC1(i1Min,1)],10.*[0 PC1(i1Min,2)],'m')
% 
% plot(10.*[0 PC2(i2Max,1)],10.*[0 PC2(i2Max,2)],'m')
% plot(10.*[0 PC2(i2Min,1)],10.*[0 PC2(i2Min,2)],'m')

idx = getOcclusionPoints2D(PC1);
plot(10.*[0 PC1(idx(1),1)],10.*[0 PC1(idx(1),2)],'m')
plot(10.*[0 PC1(idx(2),1)],10.*[0 PC1(idx(2),2)],'m')

idx = getOcclusionPoints2D(PC2);
plot(10.*[0 PC2(idx(1),1)],10.*[0 PC2(idx(1),2)],'m')
plot(10.*[0 PC2(idx(2),1)],10.*[0 PC2(idx(2),2)],'m')









