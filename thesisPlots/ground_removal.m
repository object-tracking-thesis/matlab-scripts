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
    test = pointCloud(cleanedFrames{j});
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