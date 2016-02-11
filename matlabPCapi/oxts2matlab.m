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
C = textscan(fileID,['%f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 ' ...
                    ,'%f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 ' ...
                    ,'%f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 ']); 
fclose(fileID);

pcCoord = cell2mat(C);

% lat:   latitude of the oxts-unit (deg)
% lon:   longitude of the oxts-unit (deg)
% alt:   altitude of the oxts-unit (m)
% roll:  roll angle (rad),    0 = level, positive = left side up,      range: -pi   .. +pi
% pitch: pitch angle (rad),   0 = level, positive = front down,        range: -pi/2 .. +pi/2
% yaw:   heading (rad),       0 = east,  positive = counter clockwise, range: -pi   .. +pi
% vn:    velocity towards north (m/s)
% ve:    velocity towards east (m/s)
% vf:    forward velocity, i.e. parallel to earth-surface (m/s)
% vl:    leftward velocity, i.e. parallel to earth-surface (m/s)
% vu:    upward velocity, i.e. perpendicular to earth-surface (m/s)
% ax:    acceleration in x, i.e. in direction of vehicle front (m/s^2)
% ay:    acceleration in y, i.e. in direction of vehicle left (m/s^2)
% ay:    acceleration in z, i.e. in direction of vehicle top (m/s^2)
% af:    forward acceleration (m/s^2)
% al:    leftward acceleration (m/s^2)
% au:    upward acceleration (m/s^2)
% wx:    angular rate around x (rad/s)
% wy:    angular rate around y (rad/s)
% wz:    angular rate around z (rad/s)
% wf:    angular rate around forward axis (rad/s)
% wl:    angular rate around leftward axis (rad/s)
% wu:    angular rate around upward axis (rad/s)
% pos_accuracy:  velocity accuracy (north/east in m)
% vel_accuracy:  velocity accuracy (north/east in m/s)
% navstat:       navigation status (see navstat_to_string)
% numsats:       number of satellites tracked by primary GPS receiver
% posmode:       position mode of primary GPS receiver (see gps_mode_to_string)
% velmode:       velocity mode of primary GPS receiver (see gps_mode_to_string)
% orimode:       orientation mode of primary GPS receiver (see gps_mode_to_string)