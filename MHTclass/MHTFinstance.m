% MHTFinstance class, in charge of high-level decision logic in the filter.
% Unlike the other classes, only one instance of MHTFinstance is ment to be
% used. This class is responsible for storing hypotheses, generating new
% ones, calling upon hypotheses to perform prediction, merging and pruning
% hypotheses, as well as generating new ones based on LP programming
% techniques (Murty's Method).
classdef MHTFinstance < handle
    
    properties
        hypoStorage; % This is where hypotheses will be stored
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
            
            % Setup initial hypothesis 
            
            initHypo = Hypothesis;
            initHypo.alpha = 1;
            initHypo.hypoHistory = zeros(1,scanDepth+1);
            initHypo.tracks = Tracks; % Zero tracks to being with
            
            % Create first assignment matrix (for k > 1 this will lopp instead)
%             betaFA = Model.
%             betaNT = Model.
%             assignmentMatrix =
%             
            
        end
    end
    
    %% ====================================================================
    % Internal functions
    % =====================================================================
    methods(Access = private)
        
    end
end