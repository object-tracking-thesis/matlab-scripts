function poles = poles2matlab(filepath)
% Reads file at filepath and returns cell array containing all defined
% poles
% in:
%   filepath - path to a .txt file containing the wall data
% out:
%   poles - cell array containing all defined poles

D=dir(filepath);
if isempty(D)
   error('No file found');
end

fileID = fopen(filepath);

% Read as 32-bit single precision floats
% radius, [x, y, z] (for center point), z-offset, height
C = textscan(fileID,'%f32 %f32 %f32 %f32 %f32 %f32'); 
fclose(fileID);

poles = cell2mat(C);
poles = mat2cell(poles,[ones(1,size(poles,1))],[6]);
%walls = C;