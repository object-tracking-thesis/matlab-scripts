% MHTFinstance class, in charge of high-level decision logic in the filter.
% Unlike the other classes, only one instance of MHTFinstance is ment to be
% used. This class is responsible for storing hypotheses, generating new
% ones, calling upon hypotheses to perform prediction, merging and pruning
% hypotheses, as well as generating new ones based on assignment problem
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
        fileID;
    end
    
    %% ====================================================================
    % API functions & constructor.
    % =====================================================================
    methods(Access = public)
        % Constructor 
        function this = MHTFinstance(nrHypos, scanDepth, scan)
            this.fileID = fopen('exp.txt','w');
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
%             betaFA = Model.rho;
%             betaNT = Model.spwn;
%             x = -1e10; % -inf approximation
%             
%             FAm = diag(betaFA.*ones(length(Scan.measId)));
%             FAm(FAm == 0) = x;
%             NTm = diag(betaNT.*ones(length(Scan.measId)));
%             NTm(NTm == 0) = x;
            
%            assignmentMatrix = [FAm NTm];
            
            % Get a matrix of (at max) nrHypos-best associations (e.g. 10), that we
            % use to create new hypotheses with. Each column an association
            % for the measurements (rows in the association- & assignmentMatrix)
            
            %--------------------------------------------------------------
            % TODO - Implement assignmentAlgorithm
            %--------------------------------------------------------------
            % associationMatrix = assigmentAlgorithm(assignmentMatrix, nrHypos);
            
            % DEBUG! ONLY FOR 2 measurements AND 3 hypos
            %associationMatrix = [[0 1]' [2 1]' [1 2]'];
            
            gatingMatrix = [];
            hypoGen = OptimalHypos;
            associationMatrix = hypoGen.generateHypos(scan, gatingMatrix, initHypo.tracks, this.hypoLimit);
            
            
            [~, c] = size(associationMatrix);

            % Remove empty hypos if necessary (only if c < nrHypos)
            if (c < nrHypos)
                this.hypoStorage(c+1:end) = []; 
            end
            
            % Start generating hypotheses
            for j = 1:c 
                association = associationMatrix(:,j);
                this.hypoStorage(j) = Hypothesis(initHypo, association, scan, j);
            end
            
            this.setAlphas();            
            % Since first set of hypotheses, we do no N-scan pruning. But
            % hypotheses merging should occur,             
            %--------------------------------------------------------------
            % TODO - hypotheses merging
            %--------------------------------------------------------------            
            this.setBestHypo();
        end
        
        
        function iterate(this, scan)
            % Used for new measurements when k > 1. Works similarly to the
            % constructor, expcept it generates hypoLimit^2 hypos and then
            % filters this down to hypoLimit number of hypos. 
            
            % Make room for hypos to be generated 
            this.tempStorage(1, this.hypoLimit^2) = Hypothesis;            
            %this.tempStorage;
            % Generate indexmatrix to be used when generating new hypos
            % i.e. [1 2 3;4 5 6;7 8 9] where each rownr corresponds to 
            % parenthypo, and column entries to children 
            hypoIdx = reshape(1:this.hypoLimit^2, this.hypoLimit, [])';
            
            % run through each hypothesis 
            for h = 1:length(this.hypoStorage)
                
                this.hypoStorage(h).predictTracks; % make predictiosn
                
                confInt = 0.9999999999999999999999936; % Confidence interval for gate. 
                gatingMatrix = this.getGatingMatrix(this.hypoStorage(h).tracks, scan, confInt);                
                hypoGen = OptimalHypos;
                associationMatrix = hypoGen.generateHypos(scan, gatingMatrix, this.hypoStorage(h).tracks, this.hypoLimit);
                %associationMatrix = [[0 1]' [2 1]' [1 2]'];
                                
                [~, c] = size(associationMatrix);
                % Start generating hypotheses
                for j = 1:c
                    association = associationMatrix(:,j);
                    this.tempStorage(hypoIdx(h,j)) = Hypothesis(this.hypoStorage(h), association, scan, j);
                end
            end 
                %this.tempStorage;
                % keep the nrHypos best, based on beta. Then, set alpha.
                % (Somwhere here is where N-scan pruning & merging should come in.)
                % Sorting is based on beta values
                               
                this.sortHypos();
                this.mergeHypos(); % Uses the sortHypos() function as well 
                try 
                    this.hypoStorage = this.tempStorage(1:this.hypoLimit); % Now we have final hypos
                catch e %#ok<NASGU>
                    this.hypoStorage = this.tempStorage;
                end
                this.setAlphas(); % Set the alpha values;
                this.removePoorAlphas();
                this.setBestHypo();                

                format longg
                plt =[this.hypoStorage.alpha]';
                fprintf(this.fileID, evalc('disp(plt)'));              
                fprintf(this.fileID, '\n');                                              

        end
    end
    
    %% ====================================================================
    % Internal functions
    % =====================================================================
    methods(Access = private)
        
        function gatingMat = getGatingMatrix(this, tracks, scan, Pg)
            % Generates a gating matrix. This is used to check which
            % measurements fall inside which targets gate. The returned
            % matrix is of size
            % length(scan.measurements)-by-length(tracks.track),
            % i.e. each row corresponds to a measurement and each column to
            % a track. Elements in the matrix that are equal to 1 indicate
            % that the given measurement is inside of the gate of the
            % corresponding track.
            %
            % gatingMat = getGatingMatrix(tracks, scan, Pg)
            % 
            %   tracks - Tracks object containing tracks. 
            %   scan   - Scaan object containing measurements 
            %   Pg     - Confidence interval bound.
            if isempty(tracks.trackId)
                gatingMat = [];
            else
                nrTg = length(tracks.trackId);
                nrMe = length(scan.measId);
                gatingMat = zeros(nrMe, nrTg);
                
                for tg = 1:nrTg
                    target = tracks.track(tg);
                    gatedMeasId = this.gating(scan, target, Pg);
                    gatingMat(:,tg) = gatedMeasId';
                end
            end
        end
        
        function gatedMeasId = gating(~, scan, target, Pg)
            % Support function for use in getGatingMatrix. This checks which
            % measurements fall within the gate of the target supplied to
            % the function. 
            [r, c] = size(scan.measurements);
            threshold = chi2inv(Pg, r);
            gatedMeasId = scan.measId;
            
            Sk = Model.H*target.covariance*Model.H' + Model.R;
            predictedMeas = Model.H*target.expectedValue;
            
            for m = 1:c
                d = scan.measurements(:,m) - predictedMeas;
                cal = d'*inv(Sk)*d;
                
                if (cal < threshold)
                    gatedMeasId(m) = 1;
                else
                    gatedMeasId(m) = 0;
                end
            end
            gatedMeasId(gatedMeasId == inf) = [];
            
        end
        
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
        
        function sortHypos(this)
            % Will sort hypotheses according to their beta, descending
            % order.
            [~, sortId] = sort([this.tempStorage(:).beta], 'descend');
            this.tempStorage = this.tempStorage(:,sortId);            
        end
        
        function sortHyposH(this)
            % Will sort hypotheses according to their beta, descending
            % order.
            [~, sortId] = sort([this.hypoStorage(:).beta], 'descend');
            this.hypoStorage = this.hypoStorage(:,sortId);            
        end
        
        function removePoorAlphas(this)
            this.hypoStorage = this.hypoStorage([this.hypoStorage.alpha] > 1e-20);        
        end
        
        function mergeHypos(this)
            % Function for merging hypotheses that are similar to
            % eachother. To merge two hypotheses, they need to have: i. The
            % same amount of tracks. ii. the similarity between the tracks
            % is high. Similarity is measured by the min-distance between
            % the position estimates between the targets. The merged
            % hypothesis is calculated by weighing each component with their
            % beta. parentHistory is kept from the most proable of the two
            % hypos. 
            oldSize = 1;
            newSize = 0;
            mergedStorage = []; % Temporary storage for merged pairs 
            while (oldSize > newSize) && length(this.tempStorage) > 1 % We keep iterating while we see change in nr of hypos merged
                oldSize = length(this.tempStorage); % What is the current amount of hypos?            
                for i = 1:2:length(this.tempStorage)-1 % Do pairwise comparisons between neighbors 
                    
                    left = this.tempStorage(i).tracks.track; % The left one of the pair 
                    right = this.tempStorage(i+1).tracks.track; % the right one of the pair 
                    
                    if (length(left) == length(right)) && length(left) > 0 %#ok<ISMT>
                        X = [left.expectedValue]; % Take all of the targets 
                        Y = [right.expectedValue];
                        % Investiage where min distance is. I points to
                        % which target in Y is similar to target in X
                        [D, I] = pdist2(X(1:2:3,:)',Y(1:2:3,:)','euclidean','Smallest',1); 
                        
                        totB = this.tempStorage(i).beta + this.tempStorage(i+1).beta;
                        if (mean(D) < Model.mergeThreshold) && (totB > 0)% This is our threshold for merging hypos                             
                            c = Tracks; % Empty tracks object 
                            %totB = this.tempStorage(i).beta + this.tempStorage(i+1).beta;
                            if totB == 0 % To avoid numerical instability
                                totB = 1;
                                disp('Zero beta')
                            end
                            leftB = this.tempStorage(i).beta;
                            rightB = this.tempStorage(i+1).beta;
                            % Merge the tracks according to I
                            disp('Merging!')
                            for k = 1:length(I)
                                expVal = (leftB/totB)*X(:,k) + (rightB/totB)*Y(:,I(k));
                                
                                Xcov = this.tempStorage(i).tracks.track(k).covariance;
                                Ycov = this.tempStorage(i+1).tracks.track(I(k)).covariance;
                                cov = (leftB/totB)*Xcov  + (rightB/totB)*Ycov;
                                
                                p = Posterior(expVal, cov);
                                c.addTrack(p);                                
                            end
                            h = this.tempStorage(i).copy(); % Make new hypo, based on most probable of the pair 
                            h.beta = h.beta + this.tempStorage(i+1).beta; %
                            h.tracks = c;
                            mergedStorage = [mergedStorage, h]; %#ok<AGROW>
                        else % If they are not similar enough we keep them as they are
                            mergedStorage = [mergedStorage, this.tempStorage(i), this.tempStorage(i+1)]; %#ok<AGROW>                            
                        end                        
                    else % If they don't have the same amount of targets we keep them as they are
                        mergedStorage = [mergedStorage, this.tempStorage(i), this.tempStorage(i+1)]; %#ok<AGROW>
                    end                    
                end
                
                this.tempStorage = mergedStorage; % Overwrite tempStorage with the new merged hypos
                this.sortHypos(); % Sort them to enable next round of pairwise comparisons 
                newSize = length(this.tempStorage); % Check the new size, i.e. have we managed to merge any hypos? 
                mergedStorage = [];
            end            
        end        
    end        
    
end

