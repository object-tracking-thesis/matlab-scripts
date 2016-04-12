function idx = getOcclusionPoints2D(pointCloud)
    % The functions returns the index of the two points in a pointcloud 
    % that span the occlusion-cone with respect to origin. 
    % pointCloud is a N-by-2 matrix, where each row is a data point and the
    % column values are x and y coordinates respectively. idx is a 1x2
    % matrix where each element is an index for a row in pointCloud.
    % 
    % idx = getOcclusionPoints2D(pointCloud)
    %              

    % Get the geometric center point of the pointCloud
    m = mean(pointCloud);
    % Define perpendicular line
    perpLine = [1 -m(1)/m(2)]./norm([1 -m(1)/m(2)],2);

    % project points onto perp line, using the formula proj_s-on-v =
    % s*dot(v,s)/dot(s,s)
        
    [r,c ] = size(pointCloud);
    
    top = sum(pointCloud .* repmat(perpLine,r,1),2);
    bot = sum(perpLine.*perpLine);
    
    scalingFactor = top./bot;
    projPoints = repmat(scalingFactor,1,c) .* repmat(perpLine,r,1);
    
    % Get angle for PC1 & PC2
    x = [0 1];
    phi = acos(sum(perpLine.*x)/(norm(x,2)*norm(perpLine,2)));    
    % Rotate PC1 points by phi
    R = [cos(phi) -sin(phi); sin(phi) cos(phi)];    

    projPointsRotated = projPoints*R';
    
    % From the rotated points, get the index for max and min value 
    [~, iMax] = max(projPointsRotated(:,2));
    [~, iMin] = min(projPointsRotated(:,2));
    
    idx = [iMax, iMin];
end