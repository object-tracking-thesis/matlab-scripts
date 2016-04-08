function idx = getOcclusionPoints2D(pointCloud)
    % pointCloud = [x1, y1;...
    %               x2, y2;...]

    m = mean(pointCloud);
    % Define perpendicular line
    perpLine = [1 -m(1)/m(2)]./norm([1 -m(1)/m(2)],2);

    % project points onto perp line
    projPoints = zeros(size(pointCloud));

    for k = 1:size(projPoints,1)
        projPoints(k,:) = sum(pointCloud(k,:).*perpLine) / sum(perpLine.*perpLine) .* perpLine;
    end

    % Get angle for PC1 & PC2
    x = [0 1];
    phi = acos(sum(perpLine.*x)/(norm(x,2)*norm(perpLine,2)));    
    % Rotate PC1 points by phi
    R = [cos(phi) -sin(phi); sin(phi) cos(phi)];
    

    projPointsRotated = projPoints*R';
    
    [~, iMax] = max(projPointsRotated(:,2));
    [~, iMin] = min(projPointsRotated(:,2));
    
    idx = [iMax, iMin];
end