function walls = walls2matlab(filepath)
% Reads file at filepath and returns cell array containing all defined
% walls
% in:
%   filepath - path to a .txt file containing the wall data
% out:
%   walls - cell array containing all defined walls

D=dir(filepath);
if isempty(D)
   error('No file found');
end

fileID = fopen(filepath);

% Read as 32-bit single precision floats
% [x, y, z] (for center point), [lx, ly, lz] (length), [roll, pitch, yaw]
C = textscan(fileID,'%f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32'); 
fclose(fileID);

walls = cell2mat(C);
walls = mat2cell(walls,[ones(1,size(walls,1))],[9]);
%walls = C;