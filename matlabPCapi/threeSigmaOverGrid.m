function [x1, x2, x3] = threeSigmaOverGrid(mu,P) 
    n = 100; % Number of grid points
    phi = linspace(0,2*pi,100); % Create grid in interval [0,2*pi] 4
  
    % 3 Sigma ellipse position
    x1 = repmat(mu,1,n)+1*sqrtm(P)*[cos(phi);sin(phi)];
    x2 = repmat(mu,1,n)+2*sqrtm(P)*[cos(phi);sin(phi)];
    x3 = repmat(mu,1,n)+3*sqrtm(P)*[cos(phi);sin(phi)];
end