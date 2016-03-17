%% loadGPSFile
% loads a gps.txt file and 

function [affineStart, velStart, accStart, affineEnd, velEnd, accEnd, offset] = loadGPSDir(path)
D = dir([path, '*.txt']);
Num = length(D(not([D.isdir])));
affineStart = cell(1,Num);
velStart = cell(1,Num);
accStart = cell(1,Num);
affineEnd = cell(1,Num);
velEnd = cell(1,Num);
accEnd = cell(1,Num);
offset = cell(1,Num);

    for j = 1:Num
        strcat(path,D(j).name);
        gpsData = gps2matlab(strcat(path,D(j).name));
        %get the rotation part as [1 4 7; 2 5 8; 3 6 9] of Startposition
        rotationStart = vec2mat(gpsData(:,1:9),3); 
        %get the translation part of Startposition
        translationStart = gpsData(:,10:12)';
        velStart = gpsData(:,13:15)';
        accStart = gpsData(:,16:18)';
        %get the rotation part as [1 4 7; 2 5 8; 3 6 9]
        rotationEnd = vec2mat(gpsData(:,19:27),3); 
        %get the translation part
        translationEnd = gpsData(:,28:30)';
        velEnd = gpsData(:,31:33)';
        accEnd = gpsData(:,34:36)';
        %get the offset part
        offset{j} = gpsData(:,37:39)';
        
        %add to affine cells as affine matrices
        affineStart{j} = [rotationStart translationStart; zeros(1,3) 1];
        affineEnd{j} = [rotationEnd translationEnd; zeros(1,3) 1];
    end

end