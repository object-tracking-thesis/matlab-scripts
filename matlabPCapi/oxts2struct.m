function S = oxts2struct(cellData, pathToDataformat)
% Takes cell output from oxts2matlab and converts it to a struct, specified
% by the file dataformat.txt 
% Inputs to the function are
%   cellData            - output from oxts2matlab 
%   pathToDataformat    - path to the file dataformat.txt 
   
    D=dir(pathToDataformat);
    if isempty(D)
        error('No dataformat.txt found');
    end

    fileID = fopen(pathToDataformat);
    % C is cell containing the name of the different entries 
    cellNames = textscan(fileID,'%s %*[^\n]','Delimiter',':'); 
    fclose(fileID);
    S = struct;
    for j = 1:length(cellNames{1})
        S = setfield(S,cellNames{1}{j},[]);
    end
    
    for k = 1:length(cellData)
        for j = 1:length(cellNames{1})
            S = setfield(S,cellNames{1}{j},[getfield(S,cellNames{1}{j}) cellData{k}(j)]);
        end
    end

end


