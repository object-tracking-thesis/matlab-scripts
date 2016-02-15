function translatedFrame = translateFrame(pcdXYZ, gpsOld, gpsNew, altOld, altNew)
% translates a static frame in XYZ in relation to a live frame
% according to a difference between old/new GPS and old/new altitude
% in:
%   pcdXYZ -    Nx3 matrix, [x, y, z]
%   gpsOld -    1x2 matrix, old frame position, [lat, lon] deg
%   gpsNew -    1x2 matrix new frame position, [lat, lon] deg
%   altOld -    scalar, old altitude, meter
%   altNew -    scalar, new altitude, meter
% out:
%   translatedFrame -  Nx3 matrix, [x, y, z]

xyDiff = gpsDiff2xyDiff(gpsOld, gpsNew);
zDiff = altNew - altOld;
xyzDiff = [xyDiff zDiff] %1x3

%translate for x, y and z
%e.g. if the position difference is x=+1 then
%the static frame has to subtract that difference
%from all its coordinates
pcdXYZ(:,1) = pcdXYZ(:,1) - xyzDiff(1);
pcdXYZ(:,2) = pcdXYZ(:,2) - xyzDiff(2);
pcdXYZ(:,3) = pcdXYZ(:,3) - xyzDiff(3);
translatedFrame = pcdXYZ;

end