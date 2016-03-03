function pcCoord = gps2matlab(file,varargin)
% Reads .txt file and returns 1x30 matrix where N is the amount of
% coordinates present in the pcd file. Data is x y z intensity. 
%    file     - path to txt file

D=dir(file);
if isempty(D)
   error('No file found');
end

fileID = fopen(file);

% Read as 12 64-bit double precision floats (3x3 rot matrix, lat,lon,alt)
C = textscan(fileID,'%f64 %f64 %f64 %f64 %f64 %f64 %f64 %f64 %f64 %f64 %f64 %f64'); 
fclose(fileID);

pcCoord = cell2mat(C);