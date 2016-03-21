% Class for storing Hypotheses. Not done yet
classdef Hypothesis < handle
    properties
        beta = -1;
        alpha = 0;
        hypoHistory
        hypoNr;
        tracks
        association; 
        parentHypothesis;
    end
    
    %% ====================================================================
    % API functions & constructor. These are used to create a hypothesis 
    % from a parent hypothesis, as well as performing necessary actions 
    % when creating children hypotheses. These functions include performing 
    % prediction and calculating alpha. The API functions are ment to be 
    % called from MHTFinstance object. 
    % =====================================================================
    methods (Access = public)
        
        function this = Hypothesis(parentHypothesis, association, hypoNr, gN)
            % Create active hypothesis from parentHypothesis
            if nargin == 4
                this.hypoNr = hypoNr; % Set current hypoNr
                this.tracks = Tracks; % Empty Tracks object
                this.association = association;
                this.parentHypothesis = parentHypothesis;
                this.setHypoHistory(parentHypothesis); % Update trace for this hypothesis
                this.calcBeta(gN);
                %this.updateTracks(parentHypothesis.tracks, association, Scan);
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

    %% ====================================================================
    % Internal functions used for calculating necessary parameters when
    % initiating a hypothesis. These include performing track update with
    % assigned measurements, initiating new tracks and calculating the Beta
    % value. No prediction is performed here.
    % =====================================================================
    methods (Access = public)
        
        function updateTracks(this, scan)
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
            
            % Copy over all tracks
            this.tracks = this.parentHypothesis.tracks.copy(); 
            for meas = 1:length(this.association)                
                if ismember(this.association(meas),this.parentHypothesis.tracks.trackId)
                    % The measurement is associated with an existing track
                    idx = this.association(meas); % The index of the track that the measurement is associated with                    
                    tr = this.parentHypothesis.tracks.track(idx); % Take out the track
                    m = scan.measurements(:,meas); 
                    % this is new: 
                    this.tracks.track(idx) = this.calcPosterior(m,tr); % Calculate posterior and insert it into current object 
                    
                elseif this.association(meas) == 0

                else
                    % The measurement is designated as a New Track
                    tr = this.initiateTrack(scan.measurements(:, meas));
                    this.tracks.addTrack(tr);
                end
            end            
        end
        
        function post = calcPosterior(~, measurement, prediction)
            % Calculates the new posterior based on a prediction and a
            % measurement. Calculates new posterior based on Kalman Filter
            % equations. Uses the Model class. 
            post = Posterior; % Empty Posterior object
            mu = prediction.expectedValue;
            P = prediction.covariance;
            y = measurement;
            S = Model.H * P * Model.H' + Model.R;
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
            mu = [meas(1) 0 meas(2) 0]';
            P = diag([Model.Rm, (Model.vmax/Model.kappa)^2, Model.Rm, (Model.vmax/Model.kappa)^2]);
            post.expectedValue = mu;
            post.covariance = P;
        end
        
        function calcBeta(this, gN)
            nrD = sum(length(intersect(this.parentHypothesis.tracks.trackId, this.association))); % Existing targets detected
            nrND = sum(length(setdiff(this.parentHypothesis.tracks.trackId, this.association))); % Existing targets not detected
            g0 = ((1-Model.Pd*Model.Pg)^(nrD))*((1-Model.Pd*Model.Pg)^(nrND));            
            this.beta = gN * g0 * this.parentHypothesis.alpha;
        end        
        
        function setHypoHistory(this, parentHypothesis)
            % Updates and sets the hypothesis history current hypothesis.
            % In each timestep, each generated hypothesis gets a hypothesis
            % number (hypoNr) sequentially from 1:nrHypos. The hypothesis
            % history for each history contains the hypoNumber of their
            % parents, in each timestep. So the first element in
            % hypoHistory is the hypoNr of their parent (at k-1), the
            % second of their grandparent (at k-2). This can be used in
            % N-scan pruning to detect which hypothesis are related for any
            % given depth of the hypoHistory.
            this.hypoHistory = ones(size(parentHypothesis.hypoHistory));
            this.hypoHistory(1) = parentHypothesis.hypoNr;
            this.hypoHistory(2:end) = parentHypothesis.hypoHistory(1:end-1);
        end
    end
    
end