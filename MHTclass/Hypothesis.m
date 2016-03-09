% Class for storing Hypotheses. Not done yet
classdef Hypothesis < handle
    properties
        beta
        alpha = 0;
        hypoHistory
        tracks
    end
    
    %% ====================================================================
    % API functions & constructor. These are used to create a hypothesis 
    % from a parent hypothesis, as well as performing necessary actions 
    % when creating children hypotheses. These functions include performing 
    % prediction and calculating alpha. The API functions are ment to be 
    % called from MHTFinstance object. 
    % =====================================================================
    methods (Access = public)
        
        function this = Hypothesis(parentHypothesis, association, Scan)
            % Create active hypothesis from parentHypothesis
            if nargin == 3
                % Update trace for this hypothesis
                currentTimestep = parentHypothesis.hypoHistory(1)+1;
                this.hypoHistory = ones(parentHypothesis.hypoHistory);
                this.hypoHistory(1) = currentTimestep;
                this.hypoHistory(2:end) = parentHypothesis.hypoHistory(1:end-1);
                % Use Assignment & Update tracks (Calculate Posterior for
                % each target), i.e. prediction assumed to be done in
                % previously in parent hypo. Will also set Beta.
                this.tracks = Tracks; % Empty Tracks object
                this.updateTracks(parentHypothesis.tracks, association, Scan);
                                                                                
            elseif nargin == 0
                warning('Empty hypothesis object created (no parent). Only for debugging or initiation of MHTF');
            else
                error('Wrong Nr of arguments (0 or 3 expected)');
            end
        end
        
        function setAlpha(this,totalBeta)
            this.alpha = this.beta / totalBeta;
        end
        
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

    %% ====================================================================
    % Internal functions used for calculating necessary parameters when
    % initiating a hypothesis. These include performing track update with
    % assigned measurements, initiating new tracks and calculating the Beta
    % value. No prediction is performed here.
    % =====================================================================
    methods (Access = public)
        
        function updateTracks(this,parentTracks, association, Scan)
            % Updates all tracks based on the association of measurements.
            % These include assigning each measurement to either an
            % existing track, as a FA or as a new target. If a measurement
            % is assigned to an existing track, then we expect it to be
            % assigned to an existing trackId in parentTracks. If a
            % measurement is assigned as FA, then we expect it to be
            % assigned to 0. If it does not match any of these, then the
            % measurement is assigned to a new target.
            % NOTE! If an existing target is not assigned to a measurement
            % then there occurs NO UPDATE, i.e. we use its predicted value
            % instead. 
            % Function also calculates and sets beta for current hypo. 
            
            %updatedTracks = Tracks; % Empty Tracks Object
            gN = 1;
            for k = 1:length(association)
                
                if ismember(association(k),parentTracks.trackId)
                    % The measurement is associated with an existing track
                    idx = association(k); % The index of the track that the measurement is associated with
                    tr = parentTracks.tracks.track(idx);
                    m = Scan(k);
                    a = calcPosterior(m,tr);
                    
                    this.tracks.addTrack(a);
                    
                    gN = gN * calcGn(tr, m);
                    
                elseif association(k) == 0
                    % The measurement is designated as False Alarm
                else
                    % The measurement is designated as a New Track
                    tr = this.initiateTrack(Scan.measurements(:,association(k)));
                    this.tracks.addTrack(tr);
                    
                    % Likelihood for new track, Not 100% sure that this is
                    % how it should be.                    
                    gN = gN * mvnpdf(Scan.measurements(:,association(k)),Scan.measurements(:,association(k)), Model.R); 
                end
            end
            
            this.beta = this.calcGzero(association)*gN;
        end
        
        function post = calcPosterior(measurement, prediction)
            % Calculates the new posterior based on a prediction and a
            % measurement. Calculates new posterior based on Kalman Filter
            % equations. Uses the Model class. 
            post = Posterior; % Empty Posterior object
            mu = prediction.expectedValue;
            P = prediction.covariance;
            y = measurement;
            
            S = Model.H * P * Model.H + Model.R;
            v = y - Model.H*mu;
            K = P*Model.H'*inv(S);
            
            post.expectedValue = mu + K*v;
            post.covariance = P - K*S*K';
        end
        
        function post = initiateTrack(~,meas)
            % Initiates a new track from the given measurement using single
            % point initation techniques. These necessitate the use of the
            % Model class and its variables. 
            post = Posterior; % Empty Posterior object.
            
            mu = [meas(1) meas(2) 0 0]';
            P = diag([Model.Rm, Model.Rm, (Model.vmax/Model.kappa)^2, (Model.vmax/Model.kappa)^2]);
            post.expectedValue = mu;
            post.covariance = P;
        end
        
        function beta = calcBeta(parentHypo, association, Scan)
            gZero = calcGzero(parentHypo.tracks, association);
            gN = calcGn(parentHypo.tracks, association, Scan);
            
            beta = parentHypo.alpha*gZero*gN;
        end
        
        function gZero = calcGzero(this, association)
            % Calculates the likelihood of detecting new targets and not
            % detecting existing targets. 
            % Note that Pr{association} is omitted             
            tgD = sum(association > 0); % The number of targets detected (both new and old), i.e. targets associateed to measurements 
            tgE = length(this.tracks.trackId); % Number of existing targets in current hypothesis
            nrFA = sum(association == 0); % The number of false alarms
            
            gZero = Model.rho^(nrFA)*exp(-Model.rho*Model.V)*(1-Model.Pd)^(tgE-tgD)*Model.Pd^(tgD);
        end
        
        function gN = calcGn(predictedTrack, measurement)
            % Calculates the likelihood for the measurement given the
            % (predicted) track. 
            predMu = predictedTrack.expectedValue;
            predCov = predictedTrack.covariance; 
            
            gN = mvnpdf(measurement, Model.H*predMu, Model.H*predCov*Model.H' + Model.R);
        end
    end
    
    
    
    
    
    

    
    
end