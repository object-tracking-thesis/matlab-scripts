% Posterior class for storing target posteriors. Containes two parameters,
% the expected value and the covaraince, both public. Can be instantiated
% empty or by directly setting the expected value and covariance. Used as
% building block in Tracks class.
%
% Posterior               - Emtpy posterior
% Posterior(exp, cov)     - Set posterior
% Posterior.expectedValue - the expected value of the posterior
% Posterior.covariance    - the covariance of the posterior
%
classdef Posterior < handle
    properties
        expectedValue;
        covariance;
    end
    
    methods (Access = public)
        function this = Posterior(exp, cov)
            if nargin  == 2
                this.expectedValue = exp;
                this.covariance = cov;
            elseif nargin > 0
                error('Wrong number of input arguments (2 expected)')
            end
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