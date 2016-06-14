function giw_phd_estimate = giwPHDEstimate(state, P, dof, scale, index)
%an estimate from a giw phd filter
%parametrized by:
% state     - the kinematic state (x,y,dx,dy,ddx,ddy)
% dof       - the dof inverse Wishart parameter
% scale     - the scale matrix inverse Wishart parameter

    giw_phd_estimate = struct('state',state,'P',P,'dof',dof,'scale',scale,'index',index);
end