fig = figure; fig.Name = 'raw';
fig.Position = [0 0 1800 1000];
k = 1;
plot3(data{k}(:,1), data{k}(:,2), data{k}(:,3),'.','Color',[0.7 0.7 0.7],'MarkerSize',5); hold on; axis equal
zoom(3)

dim = [.37 .45 .4 .4];
a = annotation('rectangle',dim,'Color','red');
a.LineWidth = 2;


fig.PaperPositionMode = 'auto';
path = '~/Desktop/';
print(fig,strcat(path,sprintf('%.4i',k+20)),'-dpng','-r150','-opengl');