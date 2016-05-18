function giw_meas_comp = giwMeasComp(points)
%a gaussian inverse-wishart measurement component
%parametrized by:
% points    - matrix of the points of that cluster
% center    - vector of the center coordinates
% scatter   - the scatter matrix of the points

    n = size(points,1);
    points = points';
    center = mean(points,2);
    points_sub_center = points - repmat(center,1,n);
    s_matrix = points_sub_center*points_sub_center';
    giw_meas_comp = struct('points',points,'center',center,'scatter',s_matrix,'n',n);
end