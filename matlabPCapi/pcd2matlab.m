function pcCoord = pcd2matlab(file,varargin)
% Reads .pcd file and returns Nx4 matrix where N is the amount of
% coordinates present in the pcd file. Data is x y z intensity. 
%    file     - path to pcd file
%    varargin - standard line offset is set to 11, can be changed to custom
%               value 


if length(varargin) > 1
    error('Too many inputs')
elseif length(varargin) == 1
   fprintf('Using lineoffset: %d\n',varargin{1})
   limit = varargin{1};
else
   %fprintf('Using standard lineoffset 11\n')
   limit = 11;
end

D=dir(file);
if isempty(D)
   error('No file found');
end

fileID = fopen(file);

for k=1:limit % skip uneccessary lines
    fgets(fileID); 
end
C = textscan(fileID,'%f32 %f32 %f32 %f32 %f32 %f32 %f32'); % Read as 32-bit double precision floats
fclose(fileID);

pcCoord = cell2mat(C);
