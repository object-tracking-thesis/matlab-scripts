classdef Model
    properties (Constant) 
        % Used for updating posterior 
        A = 1; % Transition Matrix 
        H = 1; % Measurement Matrix 
        Q = 1; % Covariance for Transition
        R = 1; % Covariance for Measurement 
        
        % Used for track initiation
        vmax = 4;  % Maximal (assumed) velocity
        kappa = 3; % scaling factor 
        
        % Contains detection probability and clutter specification 
        rho = 0.1; % False Measurement Density 
        V = 1;     % Volume/Area of Measurement Space
        Pd = 1;    % Detection Probability 
        
        spwn = 0.01; % New target spawning density 
    end
end