classdef UKF < handle
   properties
       predSP % sigma points for prediction
       predW  % weights for prediction
       
       updSP  % sigma points for update
       updW   % weights for update
       
       predSt % Predicted State
       predCov % Predicted Covariance 
       
       updSt  % Updated State
       updCov % Updated Covariance 
       
       Q % Motion model covariance 
       R % Measurement covariance 
       
       nMGPS % Number of MGPS
       nObsSt % Number of observed states
       nSP % Number of sigma points 
       nSt % Number of states 
       
       
   end
   
    methods
        
        function this = UKF(Q,R, nMGPS, nObsSt, nSt)            
            % Constructor, initializing noise components for system
            % Entire class assumes that W0 = 1 - n/3, i.e. Wi = 1/6;            
            this.Q = Q;
            this.R = R;
            
            this.nMGPS  = nMGPS;  % How many h-functions to we get
            this.nObsSt = nObsSt; % How many states will we observe through h-functions 
            
            this.nSt = nSt;            
            this.nSP = 1 + 2*nSt;
                        
        end
                
        function [stPred, covPred] = predictMoments(this, f)
            stPred = zeros(2,1);
            
            for j = 1:5
               stPred = stPred + f(sigPoints(:,j))*W(j); 
            end
            
            covPred = zeros(2,2);
            
            for j = 1:5
                covPred = covPred + (f(sigPoints(:,j)) - stPred)*(f(sigPoints(:,j)) - stPred)'*W(j);
            end
            
            covPred = covPred + this.Q;

        end
        
        function [stUp, covUp] = updateMoments(this, cellOfh, assignedZ)
            
            % Calculate predicted measurement
            yPred = 0;
            
            for j = 1:5
                yPred = yPred + h(SP(:,j))*W(j);
            end
            
            % Calculate Pxy
            
            Pxy = [0;
                   0];
            
            for j = 1:5
                Pxy = Pxy + (SP(:,j) - stPred)*(h(SP(:,j)) - yPred)'*W(j);
            end
            
            % Calculate S
            
            S = 0;            
            for j = 1:5
               S = S + (h(SP(:,j)) - yPred)*(h(SP(:,j)) - yPred)'*W(j);
            end
            
            S = S + this.R;
            
            K = Pxy*inv(S);
            
            stUp  = stPred + K*(Z-yPred);
            covUp = covPred - K*S*K';
            
        end
        

        function [SP, W] = generatePoints(~, st, cov)            
           SP = zeros(2, 5);
           W  = ones(1,5)*(1/6);
           W(1) = 1 - 2/3;
           
           P5 = chol(cov,'lower'); % Cholesky sqr of matrix 
           
           SP(:,1) = st; % SP-0
           
           for j = 1:2               
               SP(:,j+1)           = st + sqrt(3).*P5(:,j); % j:th colum of P5
               SP(:,j+1+2) = st - sqrt(3).*P5(:,j);               
           end                                 
        end
    end
    
end














