%% Subtract Point Cloud with Static Map
% Subtracts the static map from a given frame. The variable 'limit' is the
% threshold for filtering points (in meters). A good value is limit = 0.1;
% An additional parameter is available (K), specifying the number of neighbors
% to be used. If not specified a value of 1 is assumed.
% 
%    filteredFrame = subPC(staticMap, Frame, limit)
%    filteredFrame = subPC(staticMap, Frame, limit,K)

function filteredFrame = subPC(staticMap, Frame, limit, varargin)

    if length(varargin) == 1
        nNeigh = varargin{1};
    elseif length(varargin) > 1
        error('Too many function inputs');
    else
        nNeigh = 1;
    end
    

    [~, dist] = knnsearch(staticMap,Frame,'k',nNeigh);
    
    F = sum(dist < limit,2);
    Frame(F>0,:) = [];
    filteredFrame = Frame; 

end
