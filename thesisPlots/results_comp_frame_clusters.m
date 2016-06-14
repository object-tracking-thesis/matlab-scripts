load thesisPlots/results_scenario_desc_frame.mat
load thesisPlots/clustered_frame.mat

frame = results_scenario_desc_frame;
pointscolor = uint8(zeros(size(frame,1),3));
pointscolor(:,1:3) = 200;
pc = pointCloud(frame);
pc.Color = pointscolor;
pcshow(pc)
hold on
pcshow(clustered_frame, 'Markersize', 10)
%axis([-30 30 -20 20 -2 2])
zoom(2)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')
%legend('raw data', 'clusters')