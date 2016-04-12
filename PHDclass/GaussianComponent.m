% Class for storing each gaussian component in the PHD filter. Smallest
% buidling block used in the filter.
%
% Properties:
%   GaussianComponent.m - the mean of the gaussian
%   GaussianComponent.P - the cov of the gaussian
%   GaussianComponent.w - the associated weight
%
% Methods:
%      [] = GaussianComponent.predict()    - Performs KF prediction on the mean & cov.
%      [] = GaussianComponent.update(z)    - Performs KF update on the mean & cov.
%   lhood = GaussianComponent.calcLhood(z) - Returnes the likelhood for measurement
%                                            z conditioned on the prediction.
%


classdef GaussianComponent < handle
    properties (Access = public)
        % Component moments & weight
        m
        P
        w
        % Filter Properties
        v
        S
        K
    end
    
    methods(Access = public)
        function predict(this)
            this.m = Model.A * this.m; % m_{k|k-1}
            this.P = Model.A * this.P * Model.A' + Model.Q; % P_{k|k-1}
            
            this.v = Model.H*this.m; % innovation/predicted measurement
            this.S = Model.H*this.P*Model.H' + Model.R; % predicted measurement covariance
        end
        
        function update(this, z)
            this.K = this.P * Model.H' * inv(this.S); % Kalman gain

            this.m = this.m + this.K*(z - this.v); % Updated/Posterior mean
            this.P = this.P - this.P*Model.H'*inv(this.S)*Model.H*this.P;% Updated/Posterior covariance
        end
        
        function lhood = calcLhood(this, z)
            lhood = mvnpdf(z,this.v, this.S); % Measurement likelihood
        end
        
        % Make a copy of a handle object.
        function new = copy(this)
            % Instantiate new object of the same class.
            new = feval(class(this));
            % Copy all non-hidden properties.
            p = properties(this);
            for i = 1:length(p)
                new.(p{i}) = this.(p{i});
            end
        end
        
    end
end



