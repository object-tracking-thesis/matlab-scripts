function xyzirgbFrame = xyzirgb2matlab(filepath)
% Reads file at filepath and returns Nx7 matrix where N is the number of
% coordinates present in the pointcloud file.
% in:
%   filepath - path to pointcloud file with xyzirgb values
% out:
%   xyzirgbFrame - pointcloud Frame matrix, Nx7, x y z intensity R G B

D=dir(filepath);
if isempty(D)
   error('No file found');
end

fileID = fopen(filepath);

% Read as 32-bit single precision floats
C = textscan(fileID,'%f32 %f32 %f32 %f32 %f32 %f32 %f32'); 
fclose(fileID);

xyzirgbFrame = cell2mat(C);
