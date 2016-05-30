load data/kitti_campus_01_186_clusters.mat
%%
figure
a = 58;
b = 20;
pcshow(pointCloud(clusters{a}{b}))
hold on
plot3(-1,-1,0,'x')
testtxt = 'sensor pos.';
text(0,0,0, testtxt)
axis([-10 30 -10 10 -3 3])
hold off
xlabel('x')
ylabel('y')
zlabel('z')