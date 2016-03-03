%% loadGPSFile
% loads a gps.txt file and 

function [affine] = loadGPSDir(path)
D = dir([path, '*.txt']);
Num = length(D(not([D.isdir])));
affine = cell(1,Num);

    for j = 1:Num
        strcat(path,D(j).name);
        gpsData = gps2matlab(strcat(path,D(j).name));
        %get the rotation part as [1 4 7; 2 5 8; 3 6 9]
        rotation = vec2mat(gpsData(:,1:9),3); 
        %get the translation part
        translation = gpsData(:,10:12)';
        %padd to affine matrix
        affine{j} = [rotation translation; zeros(1,3) 1];
    end

end