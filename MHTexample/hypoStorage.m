% Class to store a generated hypothesis and associated parameters/moments 
% when using MHT. Initialized with empy constructor. 
% 
% Properties:
%   predMu     - Predicted Mean 
%   predCov    - Predicted Covariance 
%   hypo       - The association hypothesis used in current hypothesis 
%   postMu     - Posterior Mean
%   postCov    - Posterior Covariance 
%   beta       - Beta value for current hypothesis
%   alpha      - Alpha (i.e. Association Probability)
%   parentHypo - The parent hypothesis of current hypothesis 
%   hypoNr     - The hypothesis number for current hypothesis.
%                Automatically generated as random 4 digit hex, do not
%                change on your own. 
% 
classdef hypoStorage
    properties
        predMu;
        predCov;
        hypo;
        postMu;
        postCov;
        beta;
        alpha;
        parentHypo;
        hypoNr; 
    end
    methods
        function obj = hypoStorage(~)
            obj.hypoNr = randsample('0123456789abcdef',4,true);
        end
    end
end