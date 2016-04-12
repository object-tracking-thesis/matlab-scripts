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

classdef PHDinstance < handle
    properties
        birthRFS % Will be set to private later on
        componentStorage
        N = 0;% The number of targets 
    end
    
    methods(Access = public)
        
        function this = PHDinstance(weights, means, covariances)
            % Instantiate the PHD filter by generating the birth RFS
            temp(1,length(weights)) = GaussianComponent;
            for j = 1:length(weights)
                temp(j).w = weights{j};
                temp(j).m = means{j};
                temp(j).P = covariances{j};
                
                temp(j).v = Model.H * temp(j).m;
                temp(j).S = Model.H * temp(j).P * Model.H' + Model.R;
            end
            this.birthRFS = temp;            
        end
        
        function predict(this)
            % Update old components with their probability of survival, as
            % well as perform a prediction on them
            for j = 1:length(this.componentStorage)
                this.componentStorage(j).w = this.componentStorage(j).w*Model.Pd;
                this.componentStorage(j).predict();
            end
            this.componentStorage = [this.getBirthRFS(), this.componentStorage];
            
            % Set the predicted number of targets
            this.N = this.N*Model.Ps + sum([this.getBirthRFS.w]);
        end
        
        function update(this, Z)            
            % The gaussian components that we keep from prediction phase 
            predStorage(1,length(this.componentStorage)) = GaussianComponent;            
            for j = 1:length(this.componentStorage)
               predStorage(j) = this.componentStorage(j).copy();
               predStorage(j).w = predStorage(j).w * (1-Model.Pd);
            end
            
            % The new gaussian components made using measurements and our
            % predicted components 
            a = length(this.componentStorage);
            b = size(Z,2);
            updatedStorage(1,a*b) = GaussianComponent;            
            
            idx = 1;
            for m = 1:size(Z,2) % Measurement index 
               % calculate sum of weights*likelihoods & clutter parameter
               totalWL = 0;
               for j = 1:length(this.componentStorage) % component in prediction
                   gaussComp = this.componentStorage(j);%.copy();
                   totalWL = totalWL + gaussComp.w*gaussComp.calcLhood(Z(:,m));
               end
               totalWL = Model.K + Model.Pd*totalWL;
               
               for j = 1:length(this.componentStorage) % component in prediction
                   gaussComp = this.componentStorage(j).copy(); % copy over the predicted component 
                   gaussComp.w = Model.Pd * gaussComp.w * gaussComp.calcLhood(Z(:,m))/totalWL; % Calculate the new weight 
                   gaussComp.update(Z(:,m)); % KF update the gaussian component
                   updatedStorage(1,idx) = gaussComp; % Store the component
                   idx = idx + 1;
               end
               
            end
            
            this.componentStorage = [predStorage, updatedStorage]; % Store the updated PHD
            % Do some pruning 
            oldSum = sum([this.componentStorage.w]);
            
            [~,idx] = sort([this.componentStorage.w],'descend');
            this.componentStorage = this.componentStorage(idx);
            if length(this.componentStorage) > 100
                this.componentStorage = this.componentStorage(1:100);
            end
            this.componentStorage(([this.componentStorage.w] > 0.1));
            newSum = sum([this.componentStorage.w]);
            
            for k = 1:length(this.componentStorage)
               this.componentStorage(k).w = this.componentStorage(k).w * oldSum/newSum; 
            end
            
            % Set the updated nr of targets 
            this.N = this.N * (1-Model.Pd) + sum([updatedStorage.w]);
        end
        
        function bRFS = getBirthRFS(this)
            % Returns a copy of the birthRFS. This is since each Gaussian
            % Component is handle class, and we want 'new' brith RFS at
            % each timestep
            bRFS(1,length(this.birthRFS)) = GaussianComponent;
            for j = 1:length(this.birthRFS)
                bRFS(j) = this.birthRFS(j).copy();
            end
        end
        
        function gaussComps = getBestComp(this, threshold)
            idx = find([this.componentStorage.w] > threshold);
            
            gaussComps = this.componentStorage(idx);
            
            for k = 1:length(gaussComps)
                gaussComps(k) = gaussComps(k).copy();
            end
        end
        
    end
end