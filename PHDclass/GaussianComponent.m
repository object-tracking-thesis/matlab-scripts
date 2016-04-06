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
%      [] = GaussianComponent.calcPost(z)  - Performs KF update on the mean & cov.
%   lhood = GaussianComponent.calcLhood(z) - Returnes the likelhood for measurement 
%                                            z conditioned on the prediction.
%


classdef GaussianComponent < handle
    properties
        m
        P
        w
    end
    
    methods(Access = public)
        function predict(this)
            this.m = Model.
        end
        
        function calcPost(this, z)
            % TODO
        end
        function lhood = calcLhood(this, z)
            % TODO 
        end
    end
end