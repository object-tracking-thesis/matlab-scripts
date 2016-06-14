%% loadLidarDir
% loads all .pcd files for a given directory path. Returns cell containing
% the loaded frames. Is dependant on pcd2matlab function. 

function lidarData = loadLidarDir(path)
D = dir([path, '*.pcd']);
Num = length(D(not([D.isdir])));
lidarData = cell(1,Num);

    for j = 1:Num
        j
        lidarData{j} = pcd2matlab(strcat(path,D(j).name)); 
    end

end