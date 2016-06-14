% Class representing a PHD filter (instance). Handles top level logic for
% the PHD, such as predicting & updating the intensity. The filter is used
% in the follwing way:
%
% Initialize a phd filter with birthRFS
% weights = {2};
% means = {mu};
% covariances = {P};
% 
% PHD = PHDinstance(weights, means, covariances) - initalizes the filter
%                                                  with the specified birth RFS
%  [] = PHD.predict()                            - performs a prediction on the phd
%  [] = PHD.update(Z)                            - performs an update, where Z is a list of measurements 
%  Gc = PHD.getBestComp(w);                      - Returns all components that fulfill the wieght threshold w

classdef PHDinstance3 < handle
    properties
        birthRFS % Will be set to private later on
        componentStorage
        N = 0;% The number of targets 
        K
    end
    
    methods(Access = public)
        
        function this = PHDinstance3(weights, means)
            % Instantiate the PHD filter by generating the birth RFS
            temp(1,length(weights)) = CarTarget;
            for j = 1:length(weights)                
                temp(j).init(means{j}, [], [], weights{j})
            end
            this.birthRFS = temp;
            this.K = 0;
        end
        
        function predict(this)
            % Update old components with their probability of survival, as
            % well as perform a prediction on them
            this.componentStorage = [this.getBirthRFS(), this.componentStorage];
            
            for j = 1:length(this.componentStorage)
                this.componentStorage(j).weight = this.componentStorage(j).weight*0.99;
                this.componentStorage(j).predict();
            end
            
            % Set the predicted number of targets
            this.N = this.N*0.99 + sum([this.getBirthRFS.weight]);
        end
        
        function update(this, Z)            
            % The gaussian components that we keep from prediction phase 
            predStorage(1,length(this.componentStorage)) = CarTarget;            
            for j = 1:length(this.componentStorage)
               predStorage(j) = this.componentStorage(j).copy();
               predStorage(j).weight = predStorage(j).weight * (1-0.99);
            end
            
            % The new gaussian components made using measurements and our
            % predicted components 
            a = length(this.componentStorage); % Number of components 
            b = size(Z,2); % Number of measurements 
            updatedStorage(1,a*b) = CarTarget;            
            
            idx = 1;
            for m = 1:size(Z,2) % Measurement index 
               % calculate sum of weights*likelihoods & clutter parameter
               totalWL = 0;
               for j = 1:length(this.componentStorage) % component in prediction                   
                   carTarget = this.componentStorage(j).copy();
                   totalWL = totalWL + carTarget.weight*carTarget.calcLikelihood(Z{m});
               end
               totalWL = this.K + 0.99*totalWL;
               
               for j = 1:length(this.componentStorage) % component in prediction
                   carTarget = this.componentStorage(j).copy(); % copy over the predicted component 
                   carTarget.weight = 0.99 * carTarget.weight * carTarget.calcLikelihood(Z{m})/totalWL; % Calculate the new weight 
                   carTarget.update(); % KF update the gaussian component
                   updatedStorage(1,idx) = carTarget; % Store the component
                   idx = idx + 1;
               end
               
            end
            
            this.componentStorage = [predStorage, updatedStorage]; % Store the updated PHD
            % Do some pruning 
            oldSum = sum([this.componentStorage.weight]);
            
            [~,idx] = sort([this.componentStorage.weight],'descend');
            this.componentStorage = this.componentStorage(idx);
            if length(this.componentStorage) > 2 % For campus this is set to 2
                this.componentStorage = this.componentStorage(1:2);
            end
            this.componentStorage = this.componentStorage(([this.componentStorage.weight] > 0.05)); % For some reason this line fucks stuff up
            newSum = sum([this.componentStorage.weight]);
            
            for k = 1:length(this.componentStorage)
               this.componentStorage(k).weight = this.componentStorage(k).weight * oldSum/newSum; 
            end
            
            % Set the updated nr of targets 
            this.N = this.N * (1-0.99) + sum([updatedStorage.weight]);
        end
        
        function bRFS = getBirthRFS(this)
            % Returns a copy of the birthRFS. This is since each Gaussian
            % Component is handle class, and we want 'new' brith RFS at
            % each timestep
            bRFS(1,length(this.birthRFS)) = CarTarget;
            for j = 1:length(this.birthRFS)
                bRFS(j) = this.birthRFS(j).copy();
            end
        end
        
        function carTarget = getBestRect(this, threshold)
            idx = find([this.componentStorage.weight] > threshold);
            
            carTarget = this.componentStorage(idx);
            
            for k = 1:length(carTarget)
                carTarget(k) = carTarget(k).copy();
            end
        end
        
        function mergedTarget = mergeTargets(this, t1, t2)
        
        end
        
        function bool = isMerge(this,t1,t2)
           [st1, cov1] = t1.getState();
           [st2, cov2] = t2.getState();                                 
        end
        
    end
end