%TODO go from concatening vectors to preallocating the right amount of
%space whenever possible
classdef PHDfilter < handle
    properties
        birth_rfs = [];
        gaussians = [];
        F
        Q
        H
        R
        ps = 0.99;
        pd = 0.98;
        kk = 2;
        min_survival_weight = 0.00001;
        min_merge_dist = 0.1;
        max_gaussians = 50;
        number_of_targets = 0;
        index = 1;
    end
    
    methods(Access = public)
        function this = PHDfilter()
         
        end        
        
        % set the birth RFS
        function set_birth_rfs(this, means, covariances, weights)
            this.birth_rfs = [];
            for i = 1:length(weights)
                %TODO dynamic index
                gaussian_birth = ...
                    gaussianComp(means{i}, covariances{i}, weights{i}, this.index);
                this.index = this.index+1;
                this.birth_rfs = [this.birth_rfs gaussian_birth];
            end 
        end
        
        % set model parameters
        function set_model_parameters(this, F, Q, H, R)
            this.F = F;
            this.Q = Q;
            this.H = H;
            this.R = R;
        end
        
        % run predictions, for both birth RFS and existing targets
        function predict(this)
            this.gaussians = [this.gaussians this.birth_rfs];
            l = length(this.gaussians);
            for i=1:l
                this.gaussians(i).weight = this.ps*this.gaussians(i).weight;
                this.gaussians(i).mu = this.F*this.gaussians(i).mu;
                this.gaussians(i).P = ...
                    this.Q + this.F*this.gaussians(i).P*this.F';
            end
        end
        
        %update with current measurement Z
        function update(this, Z)
            curr_gaussians = [];
            for i=1:length(Z)
                weightsum = 0;
                new_gaussians = [];
                for j=1:length(this.gaussians)    

                    N = this.H*this.gaussians(j).mu;
                    S = this.R + this.H*this.gaussians(j).P*this.H';
                    K = this.gaussians(j).P*this.H'*inv(S);
                    P = (eye(size(K,1))-K*this.H)*this.gaussians(j).P;
                    w = this.pd*this.gaussians(j).weight*mvnpdf(Z(:,i), N, S);
                    mu = this.gaussians(j).mu + K*(Z(:,i)-N);
                    ind = this.gaussians(j).index;

                    weightsum = weightsum + w;
                    new_gaussians = [new_gaussians gaussianComp(mu,P,w,ind)];
                end
                for k=1:length(new_gaussians)
                    new_gaussians(k).weight = ...
                        new_gaussians(k).weight/(this.kk+weightsum);
                end
                curr_gaussians = [curr_gaussians new_gaussians];           
            end
            this.gaussians = curr_gaussians;
            this.weight_sort_gaussians;
            this.prune;
            %this.merge;
        end
        
        function weight_sort_gaussians(this)
            weightvec = zeros(length(this.gaussians),1);
            for i=1:length(this.gaussians)
                weightvec(i) = this.gaussians(i).weight;
            end
            [~, ind] = sort(weightvec(:),'descend');               
            this.gaussians = this.gaussians(ind);
        end
        
        %prune all gaussians that have a weight below a globally defined
        %threshold
        function prune(this)
            weightvec = zeros(length(this.gaussians),1);
            to_delete = [];
            for i=1:length(this.gaussians)
                weightvec(i) = this.gaussians(i).weight;
                if this.gaussians(i).weight < this.min_survival_weight
                    to_delete = [to_delete i];
                end
            end
            this.gaussians(to_delete) = [];
            weightsum_before_prune = sum(weightvec);
            this.number_of_targets = weightsum_before_prune;
            
            
            if length(this.gaussians) > this.max_gaussians             
                this.gaussians = this.gaussians(1:this.max_gaussians);
            end
            
            this.recalculate_weights(weightsum_before_prune)
        end
        
        %merging closeby targets
        function merge(this)
            i = 1;
            while i < length(this.gaussians)
                weightsum = 0;
                mu_sum = 0;       
                ind = [];
                %check all subsequent gaussians if they are closeby
                for j=(i+1):length(this.gaussians)
                    proximity = (this.gaussians(i).mu-this.gaussians(j).mu)' ...
                        *inv(this.gaussians(i).P)...
                        *(this.gaussians(i).mu-this.gaussians(j).mu);
                    if proximity < this.min_merge_dist
                        weightsum = weightsum + this.gaussians(j).weight;
                        mu_sum = mu_sum + this.gaussians(j).weight ...
                            *this.gaussians(j).mu;
                        ind = [ind j];
                    end                        
                end
                %update mu, P and weight of gaussians(i) if there were
                %closeby gaussians that can be merged
                if ~isempty(ind)
                    mu_sum = mu_sum./weightsum;
                    P_sum = 0;
                    for j=ind
                        P_sum = P_sum + this.gaussians(j).weight*( ...
                            this.gaussians(j).P + ( ...
                            mu_sum - this.gaussians(j).mu) * (...
                            mu_sum - this.gaussians(j).mu)');
                    end       
                    P_sum = P_sum./weightsum;
                    this.gaussians(i).mu = mu_sum;
                    this.gaussians(i).P = P_sum;
                    this.gaussians(i).weight = weightsum;
                    this.gaussians(ind) = [];
                    
                    %check if there are still duplicate indices left and
                    %reassign their index if need be
                    for k=(i+1):length(this.gaussians)
                        curr_ind = this.gaussians(k).index;
                        if this.gaussians(k).index == curr_ind
                            this.index = this.index+1;
                            this.gaussians(k).index = this.index;
                        end
                    end
                end
                i = i+1;
            end
        end
        
        %recalc the weights so they sum up to the same value as before
        %pruning/merging
        function recalculate_weights(this, weightsum_before)
            weightvec = zeros(length(this.gaussians),1);
            for i=1:length(this.gaussians)
                weightvec(i) = this.gaussians(i).weight;
            end
            weightsum_after = sum(weightvec);
            weight_ratio = weightsum_before/weightsum_after;
            for j = 1:length(this.gaussians)
                this.gaussians(i).weight = ...
                    this.gaussians(i).weight*weight_ratio;
            end
        end
        
        function number_of_targets = get_number_of_targets(this)
            number_of_targets = this.number_of_targets;
        end
        
        function best_estimates = get_best_estimates(this)
            n = round(this.number_of_targets);
            best_estimates = [];
            if n > 0
                best_estimates = this.gaussians(1:n);
            end
        end
    end
end