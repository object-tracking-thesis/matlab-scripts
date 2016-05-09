% Car corner definition 
%
% ii                   iii  
%  *--------------------*
%  |                    |
%  |                    |
%  |    ----->          |
%  |                    |
%  |                    |
%  *--------------------*
%  i                   iv 
%

% Define rectangle for car 
function corner = getCarCorner(measurementCluster, predictedState)
    p1 = [0 0];
    p2 = [0 0.4];
    p3 = [1 0.4];
    p4 = [0 1];
    rect = [p1; p2; p3; p4];
    
    
    rot = @(phi) [cos(phi) -sin(phi);
                  sin(phi)  cos(phi)];

    trans = @(deltaP) repmat(deltaP, 4,1);

    phi = predictedState(4); % Predicted heading
    
    deltaP = mean(measurementCluster(:,1:2)); % Mean value of point cloud L-shape
    
    modRect = rect*rot(phi)' + trans(deltaP);

    % cC is the modified rectangle, cP is the cornerPoint
    distCorner = @(cC, cP) [(sqrt((cC(1,1) - cP(1))^2 + (cC(1,2) - cP(2))^2));...
                            (sqrt((cC(2,1) - cP(1))^2 + (cC(2,2) - cP(2))^2));...
                            (sqrt((cC(3,1) - cP(1))^2 + (cC(3,2) - cP(2))^2));...
                            (sqrt((cC(4,1) - cP(1))^2 + (cC(4,2) - cP(2))^2))];
    
    [~,~,uOp] = cornerPoint(measurementCluster);
    
    c1 = uOp(1); c2 = uOp(2);
    n1 = uOp(3); n2 = uOp(4);
    
    xc = (-n1*c1 + n2*c2);
    yc = (-n2*c1 -n1*c2);

    [~, corner] = min(distCorner(modRect, [xc yc]));

end

