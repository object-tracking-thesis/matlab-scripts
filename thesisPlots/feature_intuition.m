%load data
load feature_data_plot.mat

figure;
for i = 1:length(feature_data_plot)
    subplot(2,2,i)
    cluster = pointCloud(feature_data_plot{i}');
    pointscolor=uint8(zeros(size(feature_data_plot{i}',1),3));
    pointscolor(:,1)=ceil(rand(1)*50);
    pointscolor(:,2)=ceil(rand(1)*50);
    pointscolor(:,3)=ceil(rand(1)*50);
    cluster.Color = pointscolor;
    pcshow(cluster)
    if i == 1
        title('car')
        axis([13 17.5 -1 2 -1.5 0])
    end
    if i == 2
        title('cyclist')
        axis([11 14 3 5 -1.5 0])
    end
    if i == 3
        title('pedestrian')
        axis([15 17 6 7.5 -1.5 0])
    end
    xlabel('x')
    ylabel('y')
    zlabel('z')   
end