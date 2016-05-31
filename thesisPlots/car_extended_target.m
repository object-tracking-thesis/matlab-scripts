load data/kitti_campus_01_186_clusters.mat
%% car
figure
a = 58;
b = 20;
cluster = pointCloud(clusters{a}{b});
pointscolor=uint8(zeros(size(clusters{a}{b},1),3));
pointscolor(:,1)=ceil(rand(1)*50);
pointscolor(:,2)=ceil(rand(1)*50);
pointscolor(:,3)=ceil(rand(1)*50);
cluster.Color = pointscolor;
pcshow(cluster)
hold on
plot3(-1,-1,0,'x')
testtxt = 'sensor pos.';
text(0,0,0, testtxt)
axis([-10 30 -10 10 -3 3])
hold off
xlabel('x')
ylabel('y')
zlabel('z')

%% cyclist
figure
a = 138;
b = 23;
cluster = pointCloud(clusters{a}{b});
pointscolor=uint8(zeros(size(clusters{a}{b},1),3));
pointscolor(:,1)=ceil(rand(1)*50);
pointscolor(:,2)=ceil(rand(1)*50);
pointscolor(:,3)=ceil(rand(1)*50);
cluster.Color = pointscolor;
pcshow(cluster)
hold on
plot3(-1,-1,0,'x')
testtxt = 'sensor pos.';
text(0,0,0, testtxt)
axis([-5 25 -5 10 -3 3])
zoom(1.3)
hold off
xlabel('x')
ylabel('y')
zlabel('z')