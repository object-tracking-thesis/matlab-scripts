function translatedFrame = translateFrame(pcdXYZ, gpsOld, gpsNew, altOld, altNew)
% translates the new .pcd data in XYZ in relation to some previous .pcd frame
% according to a difference between old/new GPS and old/new altitude
% in:
%   pcdXYZ -    Nx3 matrix, [x, y, z, intensity]
%   gpsOld -    1x2 matrix, old frame position, [lat, lon] deg
%   gpsNew -    1x2 matrix new frame position, [lat, lon] deg
%   altOld -    scalar, old altitude, meter
%   altNew -    scalar, new altitude, meter
% out:
%   translatedFrame -  Nx3 matrix, [x, y, z]

xyDiff = gpsDiff2xyDiff(gpsOld, gpsNew);
zDiff = altNew - altOld;
xyzDiff = [xyDiff zDiff]; %1x3

%cut the intensity values, they are not relevant here
translatedFrame = pcdXYZ(:,1:3);
%translate for x, y and z
translatedFrame(:,1) = translatedFrame(:,1) + xyzDiff(1);
translatedFrame(:,2) = translatedFrame(:,2) + xyzDiff(2);
translatedFrame(:,3) = translatedFrame(:,3) + xyzDiff(3);

end