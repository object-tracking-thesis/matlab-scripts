classdef UKF < handle
   properties
       predSP % sigma points for prediction
       predW  % weights for prediction
       
       upSP  % sigma points for update
       upW   % weights for update
       
       predSt % Predicted State
       predCov % Predicted Covariance 
       
       upSt  % Updated State
       upCov % Updated Covariance 
       
       Q % Motion model covariance 
       R % Measurement covariance 
       
       nMGPS % Number of MGPS
       nObsSt % Number of observed states
       nSP % Number of sigma points 
       nSt % Number of states 
       
       
   end
   
    methods
        
        function this = UKF(Q,R, nMGPS, nObsSt, nSt, st0, cov0)            
            % Constructor, 
            % Entire class assumes that W0 = 1 - n/3, i.e. Wi = 1/6;            
            this.Q = Q;
            this.R = R;
            
            this.nMGPS  = nMGPS;  % How many h-functions to we get
            this.nObsSt = nObsSt; % How many states will we observe through h-functions 
            
            this.nSt = nSt;       % Nr of states in statevector             
            this.nSP = 1 + 2*nSt; % Nr of sigma points to use in filter 
            
            this.upSt  = st0;  % Inital values of state for filter to start with
            this.upCov = cov0; % Inital values of covariance for filter to start with
                        
        end
                
        function [] = predictMoments(this, f)
            
            [SP, W] = this.generatePoints(this.upSt, this.upCov);

            predictedState = zeros(this.nSt,1);
            
            for j = 1:this.nSP
                predictedState = predictedState + f(SP(:,j))*W(j);
            end
                    
            
            predictedCov = zeros(this.nSt);
            
            for j = 1:this.nSP
                predictedCov = predictedCov + (f(SP(:,j)) - predictedState)*(f(SP(:,j)) - predictedState)'*W(j);
            end
            
            predictedCov = predictedCov + this.Q;
            % Assign values to class
            this.predSP = SP;
            this.predW  = W;
            
            this.predSt  = predictedState;
            this.predCov = predictedCov;

        end
        
        function [] = updateMoments(this, hCell, assignedZ)
            
            [SP, W] = this.generatePoints(this.predSt, this.predCov);
            
            %% Calculate predicted measurement
            yPred = zeros(this.nObsSt, this.nMGPS); % We have a predicted y for each mgp, where we see the observed states 
            
            for h = 1:this.nMGPS
                for j = 1:this.nSP
                    yPred(:,h) = yPred(:,h) + hCell{h}(SP(:,j))*W(j);
                end
            end
                        
            %% Calculate Pxy
            
            % Begin with calculating h(SP)- ypred matrix 
            
            yDiff = zeros(this.nObsSt*this.nMGPS, this.nSP);
            
            for j = 1:this.nSP
               yD = [];
               for h = 1:this.nMGPS
                  yD = [yD; hCell{h}(SP(:,j))-yPred(:,h)]; 
               end
               yDiff(:,j) = yD; % Each column of yDiff should contain h(SP) - ypred for each SP
            end
            
            % Calculate Pxy
            
            Pxy = zeros(this.nSt, this.nMGPS*this.nObsSt);
            
            for j = 1:this.nSP
                Pxy = Pxy + (SP(:,j) - this.predSt)*yDiff(:,j)'*W(j);
            end
            
            %% Calculate S
            
            S = zeros(this.nMGPS*this.nObsSt);
            
            for j = 1:this.nSP
                S = S + yDiff(:,j)*yDiff(:,j)'*W(j);
            end
            
            S = S + this.R*eye(this.nMGPS*this.nObsSt);
            
            %% Calculate Parameters of interest                         
            K = Pxy*inv(S);                        
            
            this.upSt  = this.predSt  + K*(assignedZ - reshape(yPred, this.nMGPS*this.nObsSt, 1));
            this.upCov = this.predCov - K*S*K';
            % Assign values to class
            this.upSP = SP;
            this.upW  = W;
            
            
        end
        

        function [SP, W] = generatePoints(this, st, cov)            
           SP = zeros(this.nSt, this.nSP);
           W  = ones(1,this.nSP)*(1/6);
           W(1) = 1 - this.nSt/3;
           
           P5 = chol(cov,'lower'); % Cholesky sqr of matrix 
           
           SP(:,1) = st; % SP-0
           
           for j = 1:this.nSt               
               SP(:,j+1)          = st + sqrt(3).*P5(:,j); % j:th colum of P5
               SP(:,j+1+this.nSt) = st - sqrt(3).*P5(:,j);               
           end                                 
        end
    end
    
end














