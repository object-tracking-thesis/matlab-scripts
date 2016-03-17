function pcCoord = oxts2matlab(file,varargin)
% Reads .txt file and returns 1x30 matrix where N is the amount of
% coordinates present in the pcd file. Data is x y z intensity. 
%    file     - path to txt file

D=dir(file);
if isempty(D)
   error('No file found');
end

fileID = fopen(file);

% Read as 30 32-bit double precision floats
C = textscan(fileID,'%f64 %f64 %f64 %f64 %f64 %f64'); 
fclose(fileID);

pcCoord = cell2mat(C);

% lat:   latitude of the oxts-unit (deg)
% lon:   longitude of the oxts-unit (deg)
% alt:   altitude of the oxts-unit (m)
% yaw:   heading (rad),       0 = east,  positive = counter clockwise, range: -pi   .. +pi
% pitch: pitch angle (rad),   0 = level, positive = front down,        range: -pi/2 .. +pi/2
% roll:  roll angle (rad),    0 = level, positive = left side up,      range: -pi   .. +pi

