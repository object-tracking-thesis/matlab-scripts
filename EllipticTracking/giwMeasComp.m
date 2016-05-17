function giw_meas_comp = giwMeasComp(points)
%a gaussian inverse-wishart measurement component
%parametrized by:
% points    - matrix of the points of that cluster
% center    - vector of the center coordinates
% scatter   - the scatter matrix of the points

    n = size(points,1);
    center = mean(points);
    points_sub_center = points - repmat(center,n,1);
    s_matrix = sum(points_sub_center'*points_sub_center);
    giw_meas_comp = struct('points',points,'center',center,'scatter',s_matrix);
end