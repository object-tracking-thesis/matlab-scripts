%% Subtract Point Cloud with Static Map
% Subtracts the static map from a given frame. The variable 'limit' is the
% threshold for filtering points (in meters). A good value is limit = 0.1;

function filteredFrame = subPC(staticMap, Frame, limit)

    [~, dist] = knnsearch(staticMap,Frame);
    Frame(dist < limit,:) = [];
    filteredFrame = Frame; 

end
