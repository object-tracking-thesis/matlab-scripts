subplot(2,1,1)
%figure;
for j=1:1
    %test = lidarData{j};
    test = pointCloud(liveFrames{j});
    pcshow(test)
    zoom(1.2)
    pause(0.3)
end
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')
set(gca,'xtick',[-70:20:70]);         
set(gca,'ytick',[-70:20:30]);
set(gca,'ztick',[-4,0,4]);

subplot(2,1,2)
%figure;
for j=1:1
    %test = lidarData{j};
    test = pointCloud(liveFrames{j});
    pcshow(test)
    zoom(1.2)
    pause(0.3)
    hold on
    for i = 1:size(walls,1)
        plotCubes(walls{i}(1:3)',walls{i}(4),walls{i}(5),walls{i}(6),walls{i}(7:9),-2,walls{i}(3)-3)
    end
    hold off
end
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')
set(gca,'xtick',[-70:20:70]);         
set(gca,'ytick',[-70:20:30]);
set(gca,'ztick',[-4,0,4]);