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
        tempStorage = Hypothesis; % This is where hypotheses will be stored before picking the best hypoLimit ones 
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
            initHypo.hypoNr = 1; % Since only one initHypo, hypoNr = 1; 
            
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
            
            %--------------------------------------------------------------
            % TODO - Implement assignmentAlgorithm
            %--------------------------------------------------------------
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

                this.hypoStorage(j) = Hypothesis(initHypo, association, Scan, j);
            end
            
            this.setAlphas();
            
            % Since first set of hypotheses, we do no N-scan pruning. But
            % hypotheses merging should occur, 
            
            %--------------------------------------------------------------
            % TODO - hypotheses merging
            %--------------------------------------------------------------
            
            this.setBestHypo();
        end
        
        
        function iterate(this, Scan)
            % Used for new measurements when k > 1. Works similarly to the
            % constructor, expcept it generates hypoLimit^2 hypos and then
            % filters this down to hypoLimit number of hypos. 
            
            % Make room for hypos to be generated 
            this.tempStorage(1, this.hypoLimit^2) = Hypothesis;            
            
            % Generate indexmatrix to be used when generating new hypos
            % i.e. [1 2 3;4 5 6;7 8 9] where each rownr corresponds to 
            % parenthypo, and column entries to children 
            hypoIdx = reshape(1:this.hypoLimit^2, this.hypoLimit, [])';
            
            % run through each hypothesis 
            for h = 1:length(this.hypoStorage)
                
                this.hypoStorage(h).predictTracks; % make prediction 
                
                associationMatrix = [[0 1]' [2 1]' [1 2]'];
                [~, c] = size(associationMatrix);
                % Start generating hypotheses
                for j = 1:c
                    association = associationMatrix(:,j);
                    this.tempStorage(hypoIdx(h,j)) = Hypothesis(this.hypoStorage(h), association, Scan, j);
                end
            end
            
                % keep the nrHypos best, based on beta. Then, set alpha.
                % (Somwhere here is where N-scan pruning & merging should come in.)
                % Sorting is based on beta values
                [~, sortId] = sort([this.tempStorage(:).beta], 'descend');
                this.tempStorage = this.tempStorage(:,sortId);
                %
                
                % SOMETHING HERE IS FUCKED UP:
                this.hypoStorage = this.tempStorage(1:this.hypoLimit); % Now we have final hypos
                this.setAlphas(); % Set the alpha values; 
                this.setBestHypo();
            
            
        end
    end
    
    %% ====================================================================
    % Internal functions
    % =====================================================================
    methods(Access = private)
        function setAlphas(this)
            % Calculate alpha for each generated hypo stored in
            % this.hypostorage. Calculates total beta by summing the betas
            % of all hypos in this.hypoStorage, and then uses this to
            % calculate 
            
            % Calculate total beta  
            totalBeta = sum([this.hypoStorage(:).beta]);

            % Set alpha 
            for j = 1:length(this.hypoStorage)
               this.hypoStorage(j).setAlpha(totalBeta) 
            end
        end
        
        function setBestHypo(this)
            allAlphas = [this.hypoStorage(:).alpha];
            a = max(allAlphas);
            bestNr = find(allAlphas == a(1));
            
            this.bestHypo = this.hypoStorage(bestNr(1));
        end
        
    end
    
    
    
end














