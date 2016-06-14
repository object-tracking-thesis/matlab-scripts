% Class encapuslating and abstracting all classes and methods specific to
% rectangular targets.
%
%       x0 - 7x1 state vector [xk, yk, vk, phik, phiDotk, wk, lk]
%       P0 - dummy, hardcoded inside of class. Can be anything.
% clusterZ - 3xN matrix
%
% Constructor
%           this = RectTarget(x0, P0)
% Methods:
%             [] = predict()
%            lik = calcLikelihood(clusterZ) 
%             [] = update()
%  [upSt, upCov] = getState()
%
% Dependencies (Classes):
%  - RectTarget.m
%    - CarIMM.m
%      - UKF.m
%    - MGPgenerator4.m
%

classdef CarTarget < handle
    properties (Access = public)
        ukfmod
        mgpGen4
        index
        weight
        f
    end
    
    methods
        %% Constructor
        function this = CarTarget
        end
        
        function init(this, x0, ~, index, weight)
            this.index = index;
            this.weight = weight;
            % THESE NEED MORE TUNING! ALL OF THEM
            % Setup motion model parameters for imm
            
            T = 0.1; % sample time 
            % model 1 (turning)
            f1 = @(st) [st(1)+T*st(3)*cos(st(4));...
                        st(2)+T*st(3)*sin(st(4));...
                        st(3);...
                        st(4) + T*st(5);...
                        st(5);
                        st(6);
                        st(7)];
            this.f = f1;
            % model 1 (turning)
               velCov1 = 0.5^2;  % velocity covariance
            phiDotCov1 = 1^2; % turningrate covariance             
                 wCov1 = 0.01^2; % width covariance
                 lCov1 = 0.02^2; % length covariance           
                 subQ1 = diag([velCov1, phiDotCov1, wCov1, lCov1]);
            gamma1 = [0 0 1 0 0 0 0 ;
                      0 0 0 0 1 0 0 ;
                      0 0 0 0 0 1 0 ;
                      0 0 0 0 0 0 1]';
            % MOTION COVARIANCE MATRIX 
            Q1 = T*gamma1*subQ1*gamma1';                                   
            
            rCov= 0.15^2; % measurement covaraince (0.1 m)            
            % initial covariance (Good for car 1 & car 2)                 
            cov0 = [ 0.0044    0.0002    0.0208   -0.0000   -0.0008    0.0031    0.0076;
                     0.0002    0.0072    0.0054    0.0022   -0.0012    0.0058   -0.0005;
                     0.0208    0.0054    0.2942    0.0014   -0.0002    0.0235    0.0295;
                    -0.0000    0.0022    0.0014    0.0012    0.0044    0.0014   -0.0002;
                    -0.0008   -0.0012   -0.0002    0.0044    0.1587   -0.0012   -0.0014;
                     0.0031    0.0058    0.0235    0.0014   -0.0012    0.0156    0.0083;
                     0.0076   -0.0005    0.0295   -0.0002   -0.0014    0.0083    0.0303];

            st0 = x0;   % inital state
                        
            this.ukfmod = UKFmod;            
            this.ukfmod.init(Q1, rCov, st0, cov0)
                                
            N = 4; % number of additional MGPS per side (totMgps = 3+2N)            
            covMGPgate = 0.1^2;
            this.mgpGen4 = MGPgenerator4(N, covMGPgate);
            
        end
        %% API functions
        function [] = predict(this)
              this.ukfmod.predictMoments(this.f);                
        end
        
        function lik = calcLikelihood(this, clusterZ)
                            
              predSt = this.ukfmod.predSt;
              
              [gatedMgpHandles1, gatedAssignedZ1] = this.mgpGen4.generate(clusterZ, predSt);                            
              assignedZo = reshape(gatedAssignedZ1', 2*size(gatedAssignedZ1,1),1);
              
              this.ukfmod.updateMoments(gatedMgpHandles1, assignedZo);
              
              yPred0 = reshape(this.ukfmod.yPred, 2*size(this.ukfmod.yPred,2),1);              
              %this.ukfmod.yPred
              %yPred0
              
              lik = mvnpdf(assignedZo, yPred0, this.ukfmod.S);
              
              if lik < 1e-1000;
                  lik = 1e-1000;
              end
        end        
        
        function [] = update(~)
            % Update has already happened in the calcLikelihood step
            % this function exists for interface purposes 
        end
        
        function [st, cov] = getState(this)
            st = this.ukfmod.upSt;
            cov = this.ukfmod.upCov; 
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
               if strcmp(p{i},'ukfmod')
                   new.(p{i}) = this.(p{i}).copy();
               end                
            end
            
        end
        
    end
end






