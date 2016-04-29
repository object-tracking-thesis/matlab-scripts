function walls = kittiwalls2matlab(filepath)
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
% 4 points, [x y z] each
C = textscan(fileID,'%f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32'); 
fclose(fileID);

height = 5;
wallPoints = cell2mat(C);
walls = zeros(size(wallPoints,1), 9);
for i = 1:size(wallPoints,1)
    walls(i,:) = fourCorners2Cube(wallPoints(i,1:3),wallPoints(i,4:6),wallPoints(i,7:9),wallPoints(i,10:12),height);
end
walls = mat2cell(walls,[ones(1,size(walls,1))],[9]);
%walls = C;