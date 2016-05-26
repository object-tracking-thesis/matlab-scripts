function giw_meas_comp = giwMeasComp(points, type)
%a gaussian inverse-wishart measurement component
%parametrized by:
% points    - matrix of the points of that cluster
% center    - vector of the center coordinates
% scatter   - the scatter matrix of the points
% n         - number of points
% type      - type of the measurement (1=clutter, 2=car, 3=bicycle, 4=ped., 5=ped.group)

    n = size(points,1);
    points = points';
    center = mean(points,2);
    points_sub_center = points - repmat(center,1,n);
    s_matrix = points_sub_center*points_sub_center';    
    giw_meas_comp = struct('points',points,'center',center,'scatter',s_matrix,'n',n,'type',type);
end