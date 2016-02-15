function xyDiff = gpsDiff2xyDiff(gpsOld, gpsNew)
% returns the x and y diff in meters between two gps positions
% in:
%   gpsOld - gps Position to compare with, 1x2 vector, [lat lon] in deg
%   gpsNew - current gps Position, 1x2 vector, [lat lon] in deg
% out:
%   xyDiff - difference between the two positions, 1x2 vector, [x y] in m

gpsOld = deg2rad(gpsOld);
gpsNew = deg2rad(gpsNew);
gpsDiff = gpsNew - gpsOld;

radius=6378.1;

%pythagorean distance in meters
y=-gpsDiff(2)*cos((gpsOld(1)+gpsNew(1))/2)*radius*1000;
x=gpsDiff(1)*radius*1000;
%d=sqrt(x*x + y*y)

xyDiff = [x y];

end

%Haversine distance is the distance along the earth's surface and not the
%straight distance
%a=sin((gpsDiff(1))/2)^2 + cos(gpsOld(1))*cos(gpsNew(1)) * sin(gpsDiff(2)/2)^2;
%c=2*atan2(sqrt(a),sqrt(1-a));
%d2=radius*c*1000