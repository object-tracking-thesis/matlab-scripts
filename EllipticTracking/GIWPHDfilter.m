classdef GIWPHDfilter < handle
    properties
        birth_rfs = [];
        giw_comps = [];
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
        function set_model_parameters(this, F, Q, H, R, T)
            this.F = F;
            this.Q = Q;
            this.H = H;
            this.R = R;
            this.T = T;
        end
        
        % run predictions, for existing targets
        function predict(this)
            %new birth_rfs indices
            for i = 1:length(this.birth_rfs)
                this.birth_rfs(i).index = this.index;
                this.index = this.index+1;
            end
            %prediction for all previous targets
            for i = 1:length(this.giw_comps)
                this.giw_comps(i).weight = this.ps*this.giw_comps(i).weight;
                this.giw_comps(i).mu = kron(this.F,eye(this.d))*this.giw_comps(i).mu;
                this.giw_comps(i).P = this.F*this.giw_comps(i).P*this.F' + this.Q;
                temp_v = this.giw_comps(i).v;
                this.giw_comps(i).v = exp(-this.T/this.tau)*this.giw_comps(i).v;
                this.giw_comps(i).V = ((this.giw_comps(i).v-this.d-1)/(temp_v-this.d-1)) .* this.giw_comps(i).V;
            end
            this.giw_comps = [this.giw_comps this.birth_rfs];
        end
        
        %update with current measurement Z and isObject vector that
        %contains the NN predictions for each of the measurements with 1
        %being clutter, 2 car, 3 cycle, 4 pedestrian
        function update(this, meas)            
            %preallocation for all updated gaussians
            n_meas = length(meas);
            n_pred = length(this.giw_comps);
            curr_giw_comps = repmat(giwComp(0,0,0,0,0,0), 1, n_meas*n_pred);
            
            %update components
            for i = 1:n_pred
                this.giw_comps(i).K = this.giw_comps(i).P*this.H';
                this.giw_comps(i).S = this.H*this.giw_comps(i).K;                
                this.giw_comps(i).z = kron(this.H,eye(this.d))*this.giw_comps(i).mu;
                %TODO handle no detection case
                %this.giw_comps(i).weight = (1-(1-exp(-this.p_gamma)*this.pd))*this.giw_comps(i).weight;
            end
            
            %update
            counter = 0;
            for i = 1:n_meas
                %when the measurement is empty
                if isempty(meas(i).points)
                	return
                end
                weightsum = 0;
                mantissas = [];
                exponents = [];
                for j = 1:n_pred                                        
                    counter = counter+1;
                    n_points = meas(i).n;
                    S = this.giw_comps(j).S + 1/n_points;                    
                    inv_S = inv(S);
                    K = this.giw_comps(j).K * inv_S;                    
                    epsilon = meas(i).center - this.giw_comps(j).z;
                    N = inv_S*(epsilon*epsilon');                    
                    mu = this.giw_comps(j).mu + kron(K,eye(this.d))*epsilon;
                    P = this.giw_comps(j).P - K*S*K';
                    v = this.giw_comps(j).v + n_points;
                    V = this.giw_comps(j).V + N + meas(i).scatter;          
                    
                    %calculate new weight-scale as a log-likelihood
                    f1 = (-this.p_gamma*log(exp(1)) + n_points*log(this.p_gamma) + log(this.pd))...
                        -(n_points*log(this.p_beta) + (this.d/2)*(n_points*log(pi) + log(n_points) + log(S)));
                    f2 = ((this.giw_comps(j).v/2)*log(det(this.giw_comps(j).V)))...
                        -((v/2)*log(det(V)));
                    f3 = (gamma_2d_log(v/2))...
                        -(gamma_2d_log(this.giw_comps(j).v/2));
                    logw_scale = f1+f2+f3;
                    [mantissa, base10_exponent] = base10_mantissa_exponent(exp(1),logw_scale);
                    
                    %the mantissa is updated with the old weight to keep
                    %track of that information
                    mantissa = mantissa*this.giw_comps(j).weight;
                    mantissas = [mantissas mantissa];
                    exponents = [exponents base10_exponent];
                    ind = this.giw_comps(j).index;                    
                    curr_giw_comps(counter) = giwComp(mu,P,v,V,this.giw_comps(j).weight,ind);
                end
                %normalizing the weights with the weightsum for the entire measurement  
                curr_giw_comps(counter-n_pred+1:counter) = ...
                    this.normalize_weights(curr_giw_comps(counter-n_pred+1:counter), exponents, mantissas);
            end
            
            this.giw_comps = curr_giw_comps;
            this.weight_sort_components;
            this.prune;
            this.merge;
            %keep only the max_comps best components
            if length(this.giw_comps) > this.max_comps             
                this.giw_comps = this.giw_comps(1:this.max_comps);
            end            
        end
        
        function giws = normalize_weights(this, giws, exponents, mantissas)
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
            for i = 1:length(giws)                    
                giws(i).weight = norm_weights(i);                
            end
        end
        
        function weight_sort_components(this)
            weightvec = zeros(length(this.giw_comps),1);
            for i=1:length(this.giw_comps)
                weightvec(i) = this.giw_comps(i).weight;
            end
            [~, ind] = sort(weightvec(:),'descend');               
            this.giw_comps = this.giw_comps(ind);
        end
        
        %prune all gaussians that have a weight below a globally defined
        %threshold
        function prune(this)
            weightvec = zeros(length(this.giw_comps),1);
            to_delete = [];
            for i=1:length(this.giw_comps)
                weightvec(i) = this.giw_comps(i).weight;
                if this.giw_comps(i).weight < this.min_survival_weight
                    to_delete = [to_delete i];
                end
            end
            this.giw_comps(to_delete) = [];
            weightsum_before_prune = sum(weightvec);
            this.number_of_targets = weightsum_before_prune;                                   
        end
        
        %merging closeby targets
        function merge(this)
            this.weight_sort_components;
            i = 1;            
            while i < length(this.giw_comps)
                weightsum = this.giw_comps(i).weight;
                mu_sum = this.giw_comps(i).weight.*this.giw_comps(i).mu; 
                P_sum = this.giw_comps(i).weight.*this.giw_comps(i).P;
                v_sum = this.giw_comps(i).weight.*this.giw_comps(i).v;
                V_sum = this.giw_comps(i).weight.*this.giw_comps(i).V;
                ind = [];
                
                %check all subsequent components if they are closeby
                for j=(i+1):length(this.giw_comps)
                    %emphasis on spatially close targets as opposed to
                    %targets with similar velocities/accelerations
                    statej = this.giw_comps(j).mu;
                    statej(1:2) = statej(1:2).*50;
                    statei = this.giw_comps(i).mu;
                    statei(1:2) = statei(1:2).*50;
                    proximity = (statej-statei)' ...
                        *inv(kron(this.giw_comps(i).P,this.giw_comps(i).V)/(this.giw_comps(i).v+3-3*this.d-2))...
                        *(statej-statei);
                    
                    if proximity < this.min_merge_dist
                        weightsum = weightsum + this.giw_comps(j).weight;
                        mu_sum = mu_sum + this.giw_comps(j).weight ...
                            *this.giw_comps(j).mu;
                        P_sum = P_sum + this.giw_comps(j).weight ...
                            *this.giw_comps(j).P;
                        v_sum = v_sum + this.giw_comps(j).weight ...
                            *this.giw_comps(j).v;
                        V_sum = V_sum + this.giw_comps(j).weight ...
                            *this.giw_comps(j).V;
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
                    this.giw_comps(i).mu = mu_sum;
                    this.giw_comps(i).P = P_sum;
                    this.giw_comps(i).v = v_sum;
                    this.giw_comps(i).V = V_sum;
                    this.giw_comps(i).weight = weightsum;                    
                    this.giw_comps(ind) = [];                 
                end                                
                i = i+1;
            end
            
            this.weight_sort_components;
            %check if there exist duplicate indices and
            %reassign their index in that case
            for i=1:(length(this.giw_comps)-1)
                curr_ind = this.giw_comps(i).index;
                for j=(i+1):length(this.giw_comps)                    
                    if this.giw_comps(j).index == curr_ind                        
                        this.giw_comps(j).index = this.index;
                        this.index = this.index+1;
                    end
                end
                i = i+1;
            end
        end
        
        %recalc the weights so they sum up to the same value as before
        %pruning/merging
        function recalculate_weights(this, weightsum_before)
            weightvec = zeros(length(this.giw_comps),1);
            for i=1:length(this.giw_comps)
                weightvec(i) = this.giw_comps(i).weight;
            end
            weightsum_after = sum(weightvec);
            weight_ratio = weightsum_before/weightsum_after;
            for j = 1:length(this.giw_comps)
                this.giw_comps(i).weight = ...
                    this.giw_comps(i).weight*weight_ratio;
                this.giw_comps(i).weight
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
                if n > length(this.giw_comps)
                    n = length(this.giw_comps);
                end
                best_estimates = this.giw_comps(1:n);
            end
        end
        
        function gaussians = get_all_gaussians(this)
            gaussians = this.giw_comps;
        end
    end
end