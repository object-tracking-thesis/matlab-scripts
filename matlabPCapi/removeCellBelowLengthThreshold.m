function cellArray = removeCellBelowLengthThreshold(cellArray, threshold)
% returns a cell array with all cells below a certain length removed
% in:
%   cellArray - array of cells, Nxcell
% out:
%   cellArray - array of cells with cells below a length threshold
%                      removed

num = 1:length(cellArray);
A = cellfun('length',cellArray);
smaller = A < threshold;

remove = smaller.*num';
remove = remove(remove > 0);

cellArray(remove) = [];

end