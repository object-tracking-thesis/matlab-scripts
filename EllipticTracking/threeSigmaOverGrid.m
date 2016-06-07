function [x1, x2, x3] = threeSigmaOverGrid(mu,P)
%example to plot a three sigma ellipse:
%   [x1,x2,x3] = threeSigmaOverGrid(mu(1:2),cov(1:2,1:2));                
%   plot(x3(1,:),x3(2,:))

    n = 100; % Number of grid points
    phi = linspace(0,2*pi,100); % Create grid in interval [0,2*pi] 4
  
    % ellipse points
    x1 = repmat(mu,1,n)+1*sqrtm(P)*[cos(phi);sin(phi)];
    x2 = repmat(mu,1,n)+2*sqrtm(P)*[cos(phi);sin(phi)];
    x3 = repmat(mu,1,n)+3*sqrtm(P)*[cos(phi);sin(phi)];
end