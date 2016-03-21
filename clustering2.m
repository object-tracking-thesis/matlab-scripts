function pcdClusters = clustering2(pcdXYZ, d, min)
% returns clusters of points that can be reached via length d edges
% in:
%   pcdXYZ  - point cloud data in XYZ, nx3 matrix
%   d       - max. distance to be considered a relevant neighbor, scalar
%   min     - minimum number of points to be a cluster in the end
% out:
%   pcdClusters - clusters of point cloud data in XYZ, cell of vectors
%

pcdXYZWork = pcdXYZ;
clusters = cell(1,250);

i = 1;
while i < (length(pcdXYZWork)+1)
    queue = pcdXYZWork(i,1:3);
    pcdXYZWork(i,:) = [];
    j=1;
    while j < (size(queue,1)+1)
        neighbors = rangesearch(pcdXYZWork, queue(j,:), d, 'NSMethod', 'kdtree');
%         if (size(neighbors{1},2) > 10)
%             size(neighbors{1},2)
%         end
        if (size(neighbors{1},2) < 10)
            pcdXYZWork(neighbors{1}',:) = [];
            j=j+1;
            continue;
        end
        if ~isempty(neighbors{1})
            ind = size(queue,1)+1:size(neighbors{1},2)+size(queue,1);
            queue(ind,:) = pcdXYZWork(neighbors{1}',1:3);
            pcdXYZWork(neighbors{1}',:) = [];
        end
        j=j+1;
    end
    clusters{i} = queue;
    i = i+1;
end

%remove all combined clusters with less than min points
clusters = removeCellBelowLengthThreshold(clusters, min);

pcdClusters = clusters;

length(clusters)

%go from indices to actual points
% pcdClusters = cell(1,length(clusters));
% for i = 1:length(clusters)
%     pcdClusters{i} = pcdXYZ(clusters{i},1:3);
% end

end