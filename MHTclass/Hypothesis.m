% Class for storing Hypotheses. Not done yet
classdef Hypothesis < handle
    properties
        beta
        alpha
        hypoHistory
        hypoMatrix
        tracks
    end
    
    methods (Access = public)
        % Create active hypothesis from parentHypothesis
        function this = hypothesis(parentHypothesis, assignment, Scan)
            if length(varargin) == 3
                % Update trace for this hypothesis
                currentTimestep = parentHypothesis.hypoHistory(1)+1;
                this.hypoHistory(1) = currentTimestep;
                this.hypoHistory(2:end) = parentHypothesis.hypoHistory(1:end-1);
                % Get tracks from parent & make prediction on them
                this.tracks = parentHypothesis.tracks;
                this.tracks = predictTracks(this.tracks);
                
                
                
                
            elseif isempty(varargin)
                
            else
                error('Wrong Nr of arguments (0 or 3 expected)');
            end
            
            
            
        end
        
        function setAlpha(this,totalBeta)
            this.alpha = this.beta / totalBeta;
        end
        
        
    end
    
    methods (Access = public)
        
        function predictTracks(this)
            % Makes bayesian prediction on tracks present in hypothesis.
            % Uses static class Model to perform update.
            for k = this.tracks.trackId;
                mu = this.tracks.track(k).expectedValue;
                P = this.tracks.track(k).covariance;
                this.tracks.track(k).expectedValue = Model.A*mu;
                this.tracks.track(k).covariance = Model.A*P*Model.A' + Model.Q;
            end
        end
    end
    
    
    
    
end