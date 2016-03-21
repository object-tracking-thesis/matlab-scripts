%% load geometrical data for static
pathWalls = 'static/walls.txt';
pathPoles = 'static/poles.txt';
pathRoadEdges = 'static/roadEdges.txt';
walls = walls2matlab(pathWalls);
poles = poles2matlab(pathPoles);
roadEdges = roadEdges2matlab(pathRoadEdges);

%% plot all walls and poles
figure;
hold on
for i = 1:length(walls)
    plotCubes(walls{i}(1:3)',walls{i}(4),walls{i}(5),walls{i}(6),walls{i}(7:9),0,0)
end
for i = 1:length(poles)
    plotCylinder(poles{i}(2:4)',poles{i}(1),poles{i}(5),poles{i}(6))
end
test = liveFrames{50};
test(:,3) = 0;
pcshow(test)