%% loadLidarDir
% loads all .txt files for a given directory path. 
% check the dataformat in oxts2matlab.m

function oxtsData = loadOxtsDir(path)
D = dir([path, '*.txt']);
Num = length(D(not([D.isdir])));
oxtsData = cell(1,Num);

    for j = 1:Num
<<<<<<< HEAD
        strcat(path,D(j).name);
=======
>>>>>>> dd22d017df26d54e7fa657586cdc7edd5c486d6c
        oxtsData{j} = oxts2matlab(strcat(path,D(j).name)); 
    end

end