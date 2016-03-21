function [points, class, varargout] = clusterLidar(lidarData, cutOff, numP)
    % Clusters lidarData into separate classes using the 'linkage' function 
    % with parameters 'single' and 'euclidean', and the 'cluster' function
    % with the parameters 'Criterion' and 'distance'. 
    %
    %   [points, class] = clusterLidar(lidarData, cutOff, numP)
    %   [points, class, Z] = clusterLidar(lidarData, cutOff, numP)
    %   [points, class, Z, c] = clusterLidar(lidarData, cutOff, numP)
    %
    % In:
    %   lidarData - Mx3 matrix, where each row is an observation
    %   cutOff    - Cutoff value to be used when clustering data, in
    %               this case the euclidean distance between clusters.
    %   numP      - Minimal numer of observations which constitute a class.
    %               Classes containing less than this value will be omitted.
    % Out:
    %   points    - The classified coordinates. Nx3 matrix. 
    %   class     - The corresponding class for each row in 'points'. Nx1
    %               matrix.
    %   Z         - The hierarchical cluster tree, to use with dendrogram
    %   c         - The nominall cluster, to use with histogram(c, length(c))
    
    Z = linkage(single(lidarData),'single','euclidean');
    c = cluster(Z,'Cutoff',cutOff,'Criterion','distance');
    % Counts the occurrences of each element in c
    y = zeros(size(c));
    for i = 1:length(c)
        y(i) = sum(c==c(i));
    end
    
    % Only keep observations which have more observations than numP
    y = y>numP;
    col = y.*c;
    lom = col == 0;
    
    lidarData(lom,:) = [];
    class = col(col > 0);
    
    points = lidarData;
    
    varargout{1} = Z;
    varargout{2} = c;

end 