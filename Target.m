% Target class, encapsulation of the underlying classes (ellipse, rectangle)
% Also in charge of initializing targets, as well as switching between the
% two types when necessary
%
%         x0 - state vector for ellipse
%         P0 - covariance vector for ellipse
%   clusterZ - 3xN matrix
%      index - indexnumber
% targetType - character indicating targetType, 'c' - car, 'e' - ellipse
%
% Constructor: 
%       this = Target(x0, P0, index)
% Methods: (this is also the order in which methods should be called)
%            [] = predict()
%           lik = calcLikelihood(clusterZ) 
%            [] = update() 
% [upSt, upCov] = getState()
%            [] = isNewType(targetType)
% 
% ===== TODO =====
% Implement missing Ellip-specific functions. Modify to fit together with
% PHD filter.

classdef Target < handle
    properties (Access = public)
        activeTarget % The active target, either EllipTarget or RectTarget
        targetType   % Target type, 'c' - car, 'e' - ellips
        index        % PHD index
        likelihood   % Value of likelihood
        weight
        NNtag
        
        upSt         % The updated state estimate
        upCov        % The updated covariance estimate   
        upv          % The updated shape dof estimate
        upV          % The updated IW scale estimate   
    end
    
    methods
        %% Constructor
        function this = Target()
            
        end
        
        function [] = init(this, x0, P0, index, weight)
            % 'e' stands for ellipse,
            %all targets in birth RFS are set as elliptical targets
            this.targetType = 'e';
            this.index = index;
            this.weight = weight;
            % TODO
            this.activeTarget = EllipTarget;            
            this.activeTarget.init(x0, P0);
            [this.upSt, this.upCov, this.upv, this.upV] = this.activeTarget.getState();
        end
        
        %% API
        function [] = isNewType(this)
            % This is the transition function, that is used to convert from
            % one target type to another. This is called AFTER UPDATE
            if strcmp(this.targetType, 'e') && strcmp(this.NNtag, 'c') % from 'e' --> 'c'
                x = this.upSt(1);
                y = this.upSt(2);
                v = sqrt(this.upSt(3)^2 + this.upSt(4)^2);
                phi = mod(atan2(this.upSt(4),this.upSt(3)),2*pi); % heading
                phiDot = 0;
                w = 1.8;
                l = 4.7;
                
                x0 = [x, y, v, phi, phiDot, w, l]'
                
                
                alpTemp = RectTarget;
                this.activeTarget = alpTemp;
                this.activeTarget
                this.activeTarget.init(x0, 0);
                this.activeTarget
                this.targetType = 'c';
            
                %TODO implement car to ellipse transition
            elseif strcmp(this.targetType, 'c') && strcmp(this.NNtag, 'e') % from 'c' --> 'e'
                x  = this.upSt(1);
                vx = this.upSt(3)*cos(this.upSt(4));
                ax = 0; % Later, when I have modified rectangular to include acceleration, these values can be set
                
                y  = this.upSt(2);
                vy = this.upSt(3)*sin(this.upSt(4));
                ay = 0; % Later, when I have modified rectangular to include acceleration, these values can be set
                
                something1 = []; % GIW moments here
                something2 = [];
                
                x0 = [x, vx, ax, y, vy, ay, something1, something2]';
                P0 = []; %... This should be set aswell
                
                this.activeTarget = EllipTarget(x0, P0);
                this.targetType = 'e';
            end
        end
        
        function [] = predict(this)
            this.activeTarget.predict();
            
            %necessary for all targets that are not detected and skip the
            %update
            [this.upSt, this.upCov, this.upv, this.upV] = this.activeTarget.getState();
        end
        
        function hasGating = hasGating(this, clusterZ)
            hasGating = this.activeTarget.gating(clusterZ);
        end
        
        function [] = updateWeightNoGating(this)
            this.weight = this.weight * (1-this.activeTarget.Pd);
        end
        
        function [mantissa, exponent] = calcLikelihood(this, clusterZ)           
            % Storage of S and yPred is delegated the internals of
            % activeTarget
            switch clusterZ.type
                case 2
                    this.NNtag = 'c';
                otherwise
                    this.NNtag = 'e';
            end
            [mantissa, exponent] = this.activeTarget.calcLikelihood(clusterZ);
        end
        
        function [] = update(this)
            this.activeTarget.update();
            [this.upSt, this.upCov, this.upv, this.upV] = this.activeTarget.getState();
        end
        
        function [upSt, upCov, upv, upV] = getState(this)
            [upSt, upCov, upv, upV] = this.activeTarget.getState();
        end
        
        function [] = setState(this, x, P, v, V)
            this.upSt = x;
            this.upCov = P;
            this.upv = v;
            this.upV = V;
        end        
        
        function new = copy(this)
            % Instantiate new object of the same class.
            new = feval(class(this));
            % Copy all non-hidden properties.
            p = properties(this);
            for i = 1:length(p)
                new.(p{i}) = this.(p{i});
            end
            
            for i = 1:length(p)
               if strcmp(p{i},'activeTarget')
                   new.(p{i}) = this.(p{i}).copy();
               end                
            end
            
        end
        
        
    end
end