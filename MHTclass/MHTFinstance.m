% MHTFinstance class, in charge of high-level decision logic in the filter.
% Unlike the other classes, only one instance of MHTFinstance is ment to be
% used. This class is responsible for storing hypotheses, generating new
% ones, calling upon hypotheses to perform prediction, merging and pruning
% hypotheses, as well as generating new ones based on LP programming
% techniques (Murty's Method).
%
% MHTFinstance(nrHypos, scanDepth, Scan) 
%
classdef MHTFinstance < handle
    
    properties
        hypoStorage = Hypothesis; % This is where hypotheses will be stored
        scanDepth; % The scan depth that will be used for purging later on
        hypoLimit; % The amount of hypotheses that will be stored for each k
        bestHypo; % The best (i.e. most probable) hypothesis in hypothesisStorage
    end
    
    %% ====================================================================
    % API functions & constructor.
    % =====================================================================
    methods(Access = public)
        
        function this = MHTFinstance(nrHypos, scanDepth, Scan)
            % Setup for the MHTF at k = 1
            this.hypoLimit = nrHypos;
            this.scanDepth = scanDepth;
            this.hypoStorage(1,nrHypos) = Hypothesis; % Preallocate hypoStorage space 
            
            % Setup initial hypothesis
            
            initHypo = Hypothesis;
            initHypo.alpha = 1;
            initHypo.hypoHistory = zeros(1,scanDepth+1);
            initHypo.tracks = Tracks; % Zero tracks to being with
            
            % Create first assignment matrix (for k > 1 this will loop for each hypothesis instead)
            betaFA = Model.rho;
            betaNT = Model.spwn;
            x = -1e10; % -inf approximation
            
            FAm = diag(betaFA.*ones(length(Scan.measId)));
            FAm(FAm == 0) = x;
            NTm = diag(betaNT.*ones(length(Scan.measId)));
            NTm(NTm == 0) = x;
            
            assignmentMatrix = [FAm NTm];
            
            % Get a matrix of (at max) nrHypos-best associations (e.g. 10), that we
            % use to create new hypotheses with. Each column an association
            % for the measurements (rows in the association- & assignmentMatrix)
            
            % TODO - Implement assignmentAlgorithm
            % associationMatrix = assigmentAlgorithm(assignmentMatrix, nrHypos);
            
            % DEBUG! ONLY FOR 2 measurements AND 3 hypos
            associationMatrix = [[0 1]' [2 1]' [1 2]'];
            
            [~, c] = size(associationMatrix);

            % Remove empty hypos if necessary (only if c < nrHypos)
            if (c < nrHypos)
                this.hypoStorage(c+1:end) = []; 
            end
            
            % Start generating hypotheses
            for j = 1:c 
                association = associationMatrix(:,j);
                this.hypoStorage(j) = Hypothesis(initHypo, association, Scan);
            end           
            
            % Calculate alpha for each generated hypo 
            totalBeta = 0;
            for j = 1:c
                totalBeta = totalBeta + this.hypoStorage(j).beta;
            end
            
            for j = 1:c
               this.hypoStorage(j).setAlpha(totalBeta) 
            end
            
            % Since first set of hypotheses, we do no N-scan pruning,
            % altough hypotheses merging should occur.
            
            % TODO - hypotheses merging
            
            % Get best hypothesis
            allAlphas = [this.hypoStorage(:).alpha];
            a = max(allAlphas);
            bestNr = find(allAlphas == a(1));
            
            this.bestHypo = this.hypoStorage(bestNr(1));
            
        end
        
        function iterate(this, Scan)
            % Used for new measurements
            
            
        end
    end
    
    %% ====================================================================
    % Internal functions
    % =====================================================================
    methods(Access = private)
        
    end
end