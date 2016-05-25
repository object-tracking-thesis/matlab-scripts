classdef NNGIWETTIMMPHDfilter < handle
    properties
        targets = [];
        birth_targets = [];
        F
        Q
        H
        R
        T
        ps = 0.98;
        pd = 0.98;
        p_gamma = 250;
        p_beta = 1;
        theta = 1;
        tau = 5;
        sigma = 2;
        d = 2;
        min_survival_weight = 0.0001;
        min_merge_dist = 0.0005;
        max_comps = 50;
        number_of_targets = 0;
        index = 1;
    end
    
    methods(Access = public)
        function this = NNGIWETTIMMPHDfilterfilter()
            
        end        
        
        % set the birth RFS
        function set_birth_rfs(this, means)
            this.targets = [];
            for i = 1:length(means)
                birth_target = Target();
                birth_target.init(means{i}, 0, this.index, 0.1);
                this.index = this.index+1;
                this.birth_targets = [this.birth_targets birth_target];
            end 
        end
        
        % run predictions, for existing targets
        function predict(this)
            %new birth_rfs indices
            for i = 1:length(this.birth_targets)
                this.birth_targets(i).index = this.index;
                this.index = this.index+1;
            end
            %prediction for all previous targets
            for i = 1:length(this.targets)
                this.targets(i).weight = this.ps*this.targets(i).weight;
                this.targets(i).predict();
            end
            this.targets = [this.targets this.birth_targets];
        end
        
        %update with current measurement Z and isObject vector that
        %contains the NN predictions for each of the measurements with 1
        %being clutter, 2 car, 3 cycle, 4 pedestrian, 5 pedestrian groups
        function update(this, meas)            
            %preallocation for all updated targets
            n_meas = length(meas);
            n_pred = length(this.targets);
            
            %TODO preallocate here
            new_targets = [];
            for i = 1:n_meas
                for j = 1:n_pred
                    target_copy = this.targets(j).copy;
                    target_copy.deepCopy;
                    new_targets = [new_targets target_copy];
                end
            end            
                   
            %likelihood update
            counter = 0;
            for i = 1:n_meas                
                mantissas = [];
                exponents = [];
                for j = 1:n_pred                                        
                    [mantissa, base10_exponent] = new_targets(j+((i-1)*n_pred)).calcLikelihood(meas(i));
                    
                    %the mantissa is updated with the old weight to keep
                    %track of that information
                    mantissa = mantissa*this.targets(j).weight;
                    
                    %save mantissas and exponents
                    mantissas = [mantissas mantissa];
                    exponents = [exponents base10_exponent];  
                    counter = counter+1;
                end
                %normalizing the weights with the weightsum for the entire measurement  
                new_targets(counter-n_pred+1:counter) = ...
                    this.normalize_weights(new_targets(counter-n_pred+1:counter), exponents, mantissas);
            end
            
            this.targets = new_targets;
            this.weight_sort_components;
            this.prune;
            
            %run the actual update on all new targets that actually
            %survived the pruning
            for i = 1:length(this.targets)
                this.targets(i).update();
            end            
            
            this.merge;
            %keep only the max_comps best components
            if length(this.targets) > this.max_comps             
                this.targets = this.targets(1:this.max_comps);
            end            
        end
        
        function targets = normalize_weights(this, targets, exponents, mantissas)
            %some exponent preprocessing because matlab will choke on all
            %exponentials above ~50 and below ~-50
            %therefore: first normalize all towards the min exponent, then
            %normalize the max to be +50 and lower all others relatively to
            %that and finally have the lowest exponents be -10
            %now matlab can safely calculate the exponentials but we still
            %retain the relative difference between the high-weighted
            %components and the very low-weighted components are irrelevant
            %anyways
            exponents_norm = exponents-repmat(min(exponents),1,length(exponents));
            exponents_norm = exponents_norm-repmat((max(exponents_norm)-50),1,length(exponents_norm));            
            exponents_norm(exponents_norm < -10) = -10;            
            expe = exp(exponents_norm);            
            weights = mantissas.*expe;
            norm_weights = weights./sum(weights);         
            for i = 1:length(targets)                    
                targets(i).weight = norm_weights(i);                
            end
        end
        
        function weight_sort_components(this)
            weightvec = zeros(length(this.targets),1);
            for i=1:length(this.targets)
                weightvec(i) = this.targets(i).weight;
            end
            [~, ind] = sort(weightvec(:),'descend');               
            this.targets = this.targets(ind);
        end
        
        %prune all gaussians that have a weight below a globally defined
        %threshold
        function prune(this)
            weightvec = zeros(length(this.targets),1);
            to_delete = [];
            for i=1:length(this.targets)
                weightvec(i) = this.targets(i).weight;
                if this.targets(i).weight < this.min_survival_weight
                    to_delete = [to_delete i];
                end
            end
            this.targets(to_delete) = [];
            weightsum_before_prune = sum(weightvec);
            this.number_of_targets = weightsum_before_prune;                                   
        end
        
        %merging closeby targets
        function merge(this)
            this.weight_sort_components;
            i = 1;            
            while i < length(this.targets)
                weightsum = this.targets(i).weight;
                [mu_sum, P_sum, v_sum, V_sum] = this.targets(i).getState;
                mu_sum = this.targets(i).weight .* mu_sum;
                P_sum = this.targets(i).weight .* P_sum;
                v_sum = this.targets(i).weight .* v_sum;
                V_sum = this.targets(i).weight .* V_sum;
                ind = [];
                
                %check all subsequent components if they are closeby
                for j=(i+1):length(this.targets)                    
                    [x_j, P_j, v_j, V_j] = this.targets(j).getState;
                    [x_i, P_i, v_i, V_i] = this.targets(i).getState;
                    proximity = (x_j-x_i)' ...
                        *inv(kron(P_i,V_i)/(v_i+3-3*this.d-2))...
                        *(x_j-x_i);
                    
                    if proximity < this.min_merge_dist
                        weightsum = weightsum + this.targets(j).weight;
                        mu_sum = mu_sum + this.targets(j).weight * x_j;
                        P_sum = P_sum + this.targets(j).weight * P_j;
                        v_sum = v_sum + this.targets(j).weight * v_j;
                        V_sum = V_sum + this.targets(j).weight * V_j;
                        ind = [ind j];
                    end                        
                end
                %update mu, P and weight of giw_comps(i) if there were
                %closeby components that could be merged
                if ~isempty(ind)                    
                    mu_sum = mu_sum./weightsum;    
                    P_sum = P_sum./weightsum;
                    v_sum = v_sum./weightsum;
                    V_sum = V_sum./weightsum;
                    this.targets(i).setState(mu_sum, P_sum, v_sum, V_sum);
                    this.targets(i).setWeight(weightsum);                  
                    this.targets(ind) = [];                 
                end                                
                i = i+1;
            end
            
            this.weight_sort_components;
            %check if there exist duplicate indices and
            %reassign their index in that case
            for i=1:(length(this.targets)-1)
                curr_ind = this.targets(i).index;
                for j=(i+1):length(this.targets)                    
                    if this.targets(j).index == curr_ind                        
                        this.targets(j).index = this.index;
                        this.index = this.index+1;
                    end
                end
                i = i+1;
            end
        end
        
        %recalc the weights so they sum up to the same value as before
        %pruning/merging
        function recalculate_weights(this, weightsum_before)
            weightvec = zeros(length(this.targets),1);
            for i=1:length(this.targets)
                weightvec(i) = this.targets(i).weight;
            end
            weightsum_after = sum(weightvec);
            weight_ratio = weightsum_before/weightsum_after;
            for j = 1:length(this.targets)
                this.targets(i).weight = ...
                    this.targets(i).weight*weight_ratio;
                this.targets(i).weight
            end
        end
        
        function number_of_targets = get_number_of_targets(this)
            number_of_targets = this.number_of_targets;
        end
        
        function best_estimates = get_best_estimates(this)            
            this.weight_sort_components;
            n = round(this.number_of_targets);
            best_estimates = [];
            if n > 0
                if n > length(this.targets)
                    n = length(this.targets);
                end
                best_estimates = this.targets(1:n);
            end
        end
        
        function gaussians = get_all_gaussians(this)
            gaussians = this.targets;
        end
    end
end