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

classdef RectTarget < handle
    properties (Access = public)
        imm
        mgpGen4
        Pd
        index
        weight
    end
    
    methods
        %% Constructor
        function this = RectTarget
        end
        
        function init(this, x0, ~, index, weight)
            this.Pd = 0.01;
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
            
            % model 2 (no turning)
            f2 = @(st) [st(1)+T*st(3)*cos(st(4));...
                        st(2)+T*st(3)*sin(st(4));...
                        st(3);...
                        st(4);...
                        st(5);
                        st(6);
                        st(7)];  
            
            % model 1 (turning)
               velCov1 = 0.5^2;  % velocity covariance
            phiDotCov1 = 0.1^2; % turningrate covariance             
                 wCov1 = 0.03^2; % width covariance
                 lCov1 = 0.05^2; % length covariance           
                 subQ1 = diag([velCov1, phiDotCov1, wCov1, lCov1]);
            gamma1 = [0 0 1 0 0 0 0 ;
                      0 0 0 0 1 0 0 ;
                      0 0 0 0 0 1 0 ;
                      0 0 0 0 0 0 1]';
            % MOTION COVARIANCE MATRIX 
            Q1 = T*gamma1*subQ1*gamma1';            
            
            % model 2 (no turning)
               velCov2 = 0.5^2;  % velocity covariance
               phiCov2 = 0.30^2; % heading covariance             
                 wCov2 = 0.03^2; % width covariance
                 lCov2 = 0.05^2; % length covariance            
                 subQ2 = diag([velCov2, phiCov2, wCov2, lCov2]);
            gamma2 = [0 0 1 0 0 0 0;
                      0 0 0 1 0 0 0;
                      0 0 0 0 0 1 0;
                      0 0 0 0 0 0 1]';            
            % MOTION COVARIANCE MATRIX
            Q2 = T*gamma2*subQ2*gamma2';
            
            rCov= 0.2^2; % measurement covaraince (0.1 m)            
            % initial covariance
            cov0 = [ 2.7646    0.0178    0.9622    0.0046   -0.0073   -0.0065    0.0155;
                     0.0178    2.4382   -0.0282    0.0169    0.0110   -0.0117    0.0013;
                     0.9622   -0.0282    0.9927   -0.0160   -0.0097    0.0019   -0.0006;
                     0.0046    0.0169   -0.0160    5.0886    2.0106   -0.0294    0.0318;
                    -0.0073    0.0110   -0.0097    2.0106    0.9950   -0.0040    0.0174;
                    -0.0065   -0.0117    0.0019   -0.0294   -0.0040    1.0110   -0.0087;
                     0.0155    0.0013   -0.0006    0.0318    0.0174   -0.0087    1.0037];
            
            % IMM parameters (mode prior)
            TPM = [0.50 0.50;
                   0.50 0.50];
            
            nObsSt = 2; % number of observed states (how many states do we measure)
            nSt = 7;    % number of states in statevector
            st0 = x0;   % inital state
                        
            this.imm = CarIMM;
            this.imm.init(nObsSt, nSt, st0, cov0, f1, f2, Q1, Q2, rCov, TPM);
                                
            N = 1; % number of additional MGPS per side (totMgps = 3+2N)
            
            covMGPgate = 0.1^2;
            this.mgpGen4 = MGPgenerator4(N, covMGPgate);
            
        end
        %% API functions
        function [] = predict(this)
              this.imm.mmPredict();
              this.Pd = 0.01;
              
              this.imm.mmUpSt = this.imm.mmPredSt;
              this.imm.mmUpCov = this.imm.mmPredCov;
  
        end
        
        function lik = calcLikelihood(this, clusterZ)
              predictedState1 = this.imm.mmPredSt(:,1);
              predictedState2 = this.imm.mmPredSt(:,2);   
              
              [gatedMgpHandles1, gatedAssignedZ1] = this.mgpGen4.generate(clusterZ.points', predictedState1);
              [gatedMgpHandles2, ~] = this.mgpGen4.generate(clusterZ.points', predictedState2);              
              
              assignedZo = reshape(gatedAssignedZ1', 2*size(gatedAssignedZ1,1),1);
              
              if isempty(gatedMgpHandles1{1}) || isempty(gatedMgpHandles2{1}) || isempty(assignedZo)
                  lik = 1;
                  
              else
                  this.imm.mmUpdate({gatedMgpHandles1, gatedMgpHandles2}, assignedZo);
                  
                  lik = this.imm.totLikelihood;
              end
        end
        
        function [] = updateWeightNoGating(this)
            this.weight = this.weight * (1-this.Pd);
        end
        
        function [] = update(~)
            % Update has already happened in the calcLikelihood step
            % this function exists for interface purposes 
        end
        
        function [st, cov, junk1, junk2] = getState(this)
            [upSt, upCov] = this.imm.getState();
             st = upSt;
            cov = upCov;
            junk1 = [];
            junk2 = [];
        end
        
        function bool = gating(this, cluster)
            predictedState1 = this.imm.mmPredSt(:,1);
            predictedState2 = this.imm.mmPredSt(:,2);
            
            p = 1.5;
            
            bool1 = carGating(cluster.points', predictedState1, p);
            bool2 = carGating(cluster.points', predictedState2, p);
            
            if bool1 || bool2
                bool = 1;
                this.Pd = 1;
            else
                bool = 0;
            end
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
               if strcmp(p{i},'imm')
                   new.(p{i}) = this.(p{i}).copy();
               end                
            end
            
        end
        
    end
end






