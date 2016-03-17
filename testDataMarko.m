%% what clusters do we want to look at?

startFrame = 80;
endFrame = 100;
filterTestClusters = cell(1,endFrame-startFrame);
for i=startFrame:endFrame
    j = i-startFrame+1;
    filterTestClusters{j} = clusters{i};
end

%% remove z axis 
for i=startFrame:endFrame
    j = i-startFrame+1;
    for k=1:length(filterTestClusters{j})
        filterTestClusters{j}{k}(:,3) = [];
    end
end

%% calculate center of each cluster
for i=startFrame:endFrame
    j = i-startFrame+1;
    for k=1:length(filterTestClusters{j})
        filterTestClusters{j}{k} = mean(filterTestClusters{j}{k});
    end
end

%% set ego position to zero and all clusters relative to that zero pos
for i=startFrame:endFrame
    j = i-startFrame+1;
    for k=1:length(filterTestClusters{j})
        filterTestClusters{j}{k} = (filterTestClusters{j}{k} - offset{i}(1:2)');
    end
end

%% convert the subcells to xy-column matrices
for i=startFrame:endFrame
    j = i-startFrame+1;
    filterTestClusters{j} = vec2mat(cell2mat(filterTestClusters{j}),2);
end

%% plot to check that everything went alright
figure
for i=1:length(filterTestClusters)
    scatter(filterTestClusters{i}(:,1), filterTestClusters{i}(:,2),'x');
    hold on
    pause(0.5)
end
hold off

%% save relevant variables to a .mat file for marko
