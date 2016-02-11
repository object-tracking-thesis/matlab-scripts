function relevantPoints = sliceSpheric(pcdSpheric, deltaAngle)
% slices the 3d point cloud horizontally and aims to detect and remove
% ground points in those slices based on the angle between consecutive
% point vectors
% variables:
%   pcdSpheric -    Nx7 matrix, 
%                   cols: r, theta(altitude), phi(azimuth),
%                   r*sin(theta), x, y, z
%   deltaAngle -   resolution for slicing horizontally (e.g. pi/720)

relevantPoints = zeros(40000,3);
n = 1;

%cycle through all angles from -pi:deltaAngle:pi
for i=-pi:deltaAngle:pi
    slice = pcdSpheric(pcdSpheric(:,3) >= i & pcdSpheric(:,3) < (i+deltaAngle),:);
    slice = sortrows(slice,4);
    %if too few points are available, skip this iteration
    if size(slice,1) < 10
        continue;
    end
    for j=1:(length(slice)-2)
        %three consecutive points
        A = [slice(j,1); slice(j,7)];
        B = [slice(j+1,1); slice(j+1,7)];
        C = [slice(j+2,1); slice(j+2,7)];
        %if there is no angle between them, their dot product should be close to 1
        if (abs(dotProductNormalized(A,B,C))+0.3) < 1
            relevantPoints(n,:) = slice(j,5:7);
            n = n+1;
        end
    end
end

%remove rows with all zeros
relevantPoints(~any(relevantPoints,2),:) = [];
size(relevantPoints)