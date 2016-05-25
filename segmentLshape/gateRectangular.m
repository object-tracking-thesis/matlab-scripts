
% [boundedCluster] = gateRectangular(cluster, cornerPoint, massVec1, massVec2, rCov)

function [boundedCluster] = gateRectangular(cluster, cP, massVec1, massVec2, gateCov)

    massVec1unit = massVec1./norm(massVec1,2);

    massVec2unit = massVec2./norm(massVec2,2);

    gatingParam = 0.5*3*sqrt(gateCov);

    adj1 = gatingParam.*massVec1unit;
    adj2 = gatingParam.*massVec2unit;

    PcLower = cP - adj1 - adj2;
    PcUpper = cP + adj1 + adj2;

    idx = ones(length(cluster),1); 
    
    for j = 1:length(cluster)
        A = ~isInRectangle(cluster(j,1:2), PcUpper, massVec1unit.*10, massVec2unit.*10);
        B =  isInRectangle(cluster(j,1:2), PcLower, massVec1unit.*10, massVec2unit.*10);
        
        idx(j) = (A && B);
    end
    
    boundedCluster = cluster(logical(idx),:);
end


function bool = isInRectangle(mz, cP, v1, v2)
    % For rectangle 1 (upper)

    AM = mz - cP;

    AB = v2;

    AD = v1;

    bool = ((0 < dot(AM,AB)) && (dot(AM,AB)< dot(AB,AB))) && ((0 < dot(AM,AD)) && (dot(AM,AD)< dot(AD,AD)));
    
end


