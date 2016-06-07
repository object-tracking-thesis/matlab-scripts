classdef EllipPHDfilter < handle
    properties
        ellipses = [];
        birth_ellipses = [];
        d = 2;
        ps = 0.98;
        min_survival_weight = 0.001;
        min_merge_dist = 0.5;
        max_comps = 50;
        number_of_targets = 0;
        index = 1;
    end
    
    methods(Access = public)
        function this = EllipPHDfilterfilter()
            
        end        
        
        % set the birth RFS
        function set_birth_rfs(this, means)
            this.ellipses = [];
            for i = 1:length(means)
                birth_target = EllipTarget();
                birth_target.init(means{i}, 0, this.index, 0.1);
                this.index = this.index+1;
                this.birth_ellipses = [this.birth_ellipses birth_target];
            end 
            this.ellipses = this.birth_ellipses(1:7);
            this.birth_ellipses(1:7) = [];
        end
        
        % run predictions, for existing targets
        function predict(this)
            %new birth_rfs indices
            for i = 1:length(this.birth_ellipses)
                this.birth_ellipses(i).index = this.index;
                this.index = this.index+1;
            end
            %prediction for all previous targets
            for i = 1:length(this.ellipses)
                this.ellipses(i).weight = this.ps*this.ellipses(i).weight;
                this.ellipses(i).predict();
            end
            this.ellipses = [this.ellipses this.birth_ellipses];
        end
        
        %update with current measurement Z and isObject vector that
        %contains the NN predictions for each of the measurements with 1
        %being clutter, 2 car, 3 cycle, 4 pedestrian, 5 pedestrian groups
        function update(this, meas)
            length(this.ellipses)
            ind = [];
            for i = 1:length(meas)
                isMeas = 0;
                for j = 1:length(this.ellipses)                    
                    isMeas = this.ellipses(j).gating(meas(i));
                    if isMeas
                        break;
                    end
                end
                %when the measurement was not in any gate              
                if ~isMeas
                    ind = [ind i];
                end
            end
            meas(ind) = [];
            
            %set pd dynamically depending on whether a measurement was
            %received inside the gate of that target and update the weight
            ellipses_no_meas = [];
            for i = 1:length(this.ellipses)
                hasMeas = 0;
                for j = 1:length(meas)
                    hasMeas = this.ellipses(i).gating(meas(j));
                    if hasMeas
                        break;
                    end
                end
                %when there was no measurement within the gate              
                if ~hasMeas
                    this.ellipses(i).updateWeightNoGating;
                    ellipse_copy = this.ellipses(i).copy;
                    ellipses_no_meas = [ellipses_no_meas ellipse_copy];                    
                end
            end
            
            %preallocation for all targets to be updated
            n_meas = length(meas);
            n_pred = length(this.ellipses);
            
            %TODO preallocate here
            new_ellipses = [];
            for i = 1:n_meas
                for j = 1:n_pred
                    ellipse_copy = this.ellipses(j).copy;
                    new_ellipses = [new_ellipses ellipse_copy];
                end
            end            
                   
            %likelihood update
            counter = 0;
            for i = 1:n_meas                
                mantissas = [];
                exponents = [];                
                for j = 1:n_pred
                    [mantissa, base10_exponent] = new_ellipses(j+((i-1)*n_pred)).calcLikelihood(meas(i));
                    
                    %the mantissa is updated with the old weight to keep
                    %track of that information
                    mantissa = mantissa*this.ellipses(j).weight;
                    
                    %save mantissas and exponents
                    mantissas = [mantissas mantissa];
                    exponents = [exponents base10_exponent];  
                    counter = counter+1;
                end
                %normalizing the weights with the weightsum for the entire measurement  
                new_ellipses(counter-n_pred+1:counter) = ...
                    this.normalize_weights(new_ellipses(counter-n_pred+1:counter), exponents, mantissas);              
            end            
            
            this.ellipses = new_ellipses;
            this.weight_sort_components;
            this.prune;            
            
            %run the actual update on all new ellipses that actually
            %survived the pruning
            for i = 1:length(this.ellipses)
                this.ellipses(i).update();
            end        
                      
            this.ellipses = [this.ellipses  ellipses_no_meas];
            this.prune;                             
            
            this.merge;
            %keep only the max_comps best components
            if length(this.ellipses) > this.max_comps             
                this.ellipses = this.ellipses(1:this.max_comps);
            end      
        end        
        
        function ellipses = normalize_weights(this, ellipses, exponents, mantissas)
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
            mantissas(mantissas == 0) = 0.1;
            weights = mantissas.*expe;
            norm_weights = weights./sum(weights);         
            for i = 1:length(ellipses)                    
                ellipses(i).weight = norm_weights(i);                
            end
        end
        
        function weight_sort_components(this)
            weightvec = zeros(length(this.ellipses),1);
            for i=1:length(this.ellipses)
                weightvec(i) = this.ellipses(i).weight;
            end
            [~, ind] = sort(weightvec(:),'descend');               
            this.ellipses = this.ellipses(ind);
        end
        
        %prune all gaussians that have a weight below a globally defined
        %threshold
        function prune(this)
            weightvec = zeros(length(this.ellipses),1);
            to_delete = [];
            for i=1:length(this.ellipses)
                weightvec(i) = this.ellipses(i).weight;
                if this.ellipses(i).weight < this.min_survival_weight | abs(this.ellipses(i).x(5)) > 5
                    to_delete = [to_delete i];
                end
            end
            this.ellipses(to_delete) = [];
            weightsum_before_prune = sum(weightvec);
            this.number_of_targets = weightsum_before_prune;                                   
        end
        
        %merging closeby targets
        function merge(this)
            this.weight_sort_components;
            i = 1;            
            while i < length(this.ellipses)
                weightsum = this.ellipses(i).weight;
                [mu_sum, P_sum, v_sum, V_sum] = this.ellipses(i).getState;
                mu_sum = this.ellipses(i).weight .* mu_sum;
                P_sum = this.ellipses(i).weight .* P_sum;
                v_sum = this.ellipses(i).weight .* v_sum;
                V_sum = this.ellipses(i).weight .* V_sum;
                ind = [];
                
                %check all subsequent components if they are closeby
                for j=(i+1):length(this.ellipses)                    
                    [x_j, P_j, v_j, V_j] = this.ellipses(j).getState;
                    [x_i, P_i, v_i, V_i] = this.ellipses(i).getState;
                    proximity = (x_j-x_i)' ...
                        *inv(kron(P_i,V_i)/(v_i+3-3*this.d-2))...
                        *(x_j-x_i);
                    
                    if proximity < this.min_merge_dist
                        weightsum = weightsum + this.ellipses(j).weight;
                        mu_sum = mu_sum + this.ellipses(j).weight * x_j;
                        P_sum = P_sum + this.ellipses(j).weight * P_j;
                        v_sum = v_sum + this.ellipses(j).weight * v_j;
                        V_sum = V_sum + this.ellipses(j).weight * V_j;
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
                    this.ellipses(i).setState(mu_sum, P_sum, v_sum, V_sum);
                    this.ellipses(i).weight = weightsum;                  
                    this.ellipses(ind) = [];                 
                end                                
                i = i+1;
            end
            
            this.weight_sort_components;
            %check if there exist duplicate indices and
            %reassign their index in that case
            for i=1:(length(this.ellipses)-1)
                curr_ind = this.ellipses(i).index;
                for j=(i+1):length(this.ellipses)                    
                    if this.ellipses(j).index == curr_ind                        
                        this.ellipses(j).index = this.index;
                        this.index = this.index+1;
                    end
                end
                i = i+1;
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
                if n > length(this.ellipses)
                    n = length(this.ellipses);
                end
                best_estimates = this.ellipses(1:n);
            end
        end
        
        function gaussians = get_all_gaussians(this)
            gaussians = this.ellipses;
        end
    end
end