function pcdSpheric = xyzToSpheric(pcdXYZ)
% converts .pcd data in XYZ to .pcd data in Spheric
%
% variables:
%   pcdXYZ -    Nx4 matrix, 
%               cols: x, y, z, intensity
%   pcdSpheric -    Nx7 matrix, 
%                   cols: r, theta(altitude), phi(azimuth),
%                   r*sin(theta), x, y, z

n = length(pcdXYZ);
pcdSpheric = zeros(n,7);

%convert to spheric coordinates
for i=1:n
    pcdSpheric(i,1) = sqrt(pcdXYZ(i,1)^2 + pcdXYZ(i,2)^2 + pcdXYZ(i,3)^2); %r
    pcdSpheric(i,2) = acos(pcdXYZ(i,3)/pcdSpheric(i,1)); %theta = altitude
    pcdSpheric(i,3) = atan2(pcdXYZ(i,2),pcdXYZ(i,1)); %phi = azimuth [-pi,pi]
    pcdSpheric(i,4) = pcdSpheric(i,1)*sin(pcdSpheric(i,2)); %r*cos(theta)
    %save the original coordinates, too
    pcdSpheric(i,5) = pcdXYZ(i,1);
    pcdSpheric(i,6) = pcdXYZ(i,2);
    pcdSpheric(i,7) = pcdXYZ(i,3);
end