function pcdClusters = clustering(pcdXYZ, d, min)
% returns clusters of points that can be reached via length d edges
% in:
%   pcdXYZ  - point cloud data in XYZ, nx3 matrix
%   d       - max. distance to be considered a relevant neighbor, scalar
%   min     - minimum number of points to be a cluster in the end
% out:
%   pcdClusters - clusters of point cloud data in XYZ, cell of vectors
%
%TODO: ugly (and ineefficient) code

clusters = rangesearch(pcdXYZ, pcdXYZ, d);

%remove all miniclusters in the beginning (makes it a lot faster)
clusters = removeCellBelowLengthThreshold(clusters, 3);

%recursively combine all subsets that have at least one point in common
i = 1;
length(clusters)
while i < (length(clusters)-1)
    j = 1;
    while j < length(clusters)-1
        %don't compare with itself
        if i==j
            j = j+1;
            continue;
        end
        %when they have no elements in common, proceed
        %or ismember
        if ~ismembc(clusters{i},clusters{j})
            j = j+1;
        else
            %combine their unique elements
            clusters{i} = union(clusters{i},clusters{j});
            clusters(j) = [];
        end
    end
    i = i+1;
end

%remove all combined clusters with less than min points
clusters = removeCellBelowLengthThreshold(clusters, min);

%go from indices to actual points
pcdClusters = cell(1,length(clusters));
for i = 1:length(clusters)
    pcdClusters{i} = pcdXYZ(clusters{i},1:3);
end

end