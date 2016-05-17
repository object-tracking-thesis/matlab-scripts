classdef GIWPHDfilter < handle
    properties
        birth_rfs = [];
        giw_comps = [];
        F
        Q
        H
        R
        ps = 1.0;
        pd = 1.0;
        p_gamma = 100;
        p_beta = 0;
        min_survival_weight = 0.00001;
        min_merge_dist = 0.5;
        max_gaussians = 50;
        number_of_targets = 0;
        index = 1;
    end
    
    methods(Access = public)
        function this = PHDfilter()
         
        end        
        
        % set the birth RFS
        function set_birth_rfs(this, means, covariances, dofs, scale_matrices, weights)
            this.birth_rfs = [];
            for i = 1:length(weights)
                gaussian_birth = ...
                    giwComp(means{i}, covariances{i}, dofs{i}, scale_matrices{i}, weights{i}, this.index);
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
%                 %test plot 3sigma for the predictions
%                 [x1,x2,x3] = threeSigmaOverGrid(this.gaussians(i).mu(1:2),this.gaussians(i).P(1:2,1:2));                
%                 plot(x3(1,:),x3(2,:),' --k')
%                 pause(1)        
%                 hold on
            end
        end
        
        %update with current measurement Z and isObject vector that
        %contains the NN predictions for each of the measurements with 1
        %being clutter, 2 car, 3 cycle, 4 pedestrian
        function update(this, Z, isObject)
            %preallocation for all updated gaussians
            n_meas = size(Z,2);
            n_pred = length(this.gaussians);
            curr_gaussians = repmat(gaussianComp(0,0,0,0), 1, n_meas*n_pred);
            
            %cardinality of the clutter set
            n_objects = sum(isObject==2 | isObject==3 | isObject==4);
            if n_objects > 0
                this.kk = (n_meas-n_objects);
            end
            
            counter = 0;
            for i=1:size(Z,2)
                weightsum = 0;
                for j=1:length(this.gaussians)    
                    counter = counter+1;
                    N = this.H*this.gaussians(j).mu;
                    S = this.R + this.H*this.gaussians(j).P*this.H';
                    K = this.gaussians(j).P*this.H'*inv(S);
                    P = (eye(size(K,1))-K*this.H)*this.gaussians(j).P;
                    w = this.pd*this.gaussians(j).weight*mvnpdf(Z(:,i), N, S);
                    if isObject(i)==2 || isObject(i)==3 || isObject(i)==4
                        w = w*100;
                    end
                    mu = this.gaussians(j).mu + K*(Z(:,i)-N);
                    ind = this.gaussians(j).index;

                    weightsum = weightsum + w;
                    curr_gaussians(counter) = gaussianComp(mu,P,w,ind);
                end     
                %adjust the weights of the current batch of updates
                for k=0:length(this.gaussians)-1
                    curr_gaussians(counter-k).weight = ...
                        curr_gaussians(counter-k).weight/(this.kk+weightsum);
                end          
            end
            
            this.gaussians = curr_gaussians;
            this.weight_sort_gaussians;
            this.prune;
            this.merge;
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
                weightsum = this.gaussians(i).weight;
                mu_sum = this.gaussians(i).weight.*this.gaussians(i).mu;       
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
                    P_sum = this.gaussians(i).weight.*this.gaussians(i).P; 
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
                end
                
                %check if there are still duplicate indices left and
                %reassign their index if need be
                for j=(i+1):length(this.gaussians)
                    curr_ind = this.gaussians(i).index;
                    if this.gaussians(j).index == curr_ind
                        this.index = this.index+1;
                        this.gaussians(j).index = this.index;
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
                if n > length(this.gaussians)
                    n = length(this.gaussians)
                end
                best_estimates = this.gaussians(1:n);
            end
        end
        
        function gaussians = get_all_gaussians(this)
            gaussians = this.gaussians;
        end
    end
end