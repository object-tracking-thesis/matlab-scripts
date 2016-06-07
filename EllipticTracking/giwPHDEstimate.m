function giw_phd_estimate = giwPHDEstimate(state, dof, scale, index)
%an estimate from a giw phd filter
%parametrized by:
% state     - the kinematic state (x,y,dx,dy,ddx,ddy)
% dof       - the dof inverse Wishart parameter
% scale     - the scale matrix inverse Wishart parameter

    giw_phd_estimate = struct('state',state,'dof',dof,'scale',scale,'index',index);
end