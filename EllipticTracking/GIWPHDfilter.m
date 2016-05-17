classdef GIWPHDfilter < handle
    properties
        birth_rfs = [];
        giw_comps = [];
        F
        Q
        H
        R
        T
        ps = 1.0;
        pd = 1.0;
        p_gamma = 100;
        p_beta = 0;
        theta = 1;
        tau = 5;
        sigma = 2;
        d = 2;
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
        function set_model_parameters(this, F, Q, H, R, T)
            this.F = F;
            this.Q = Q;
            this.H = H;
            this.R = R;
            this.T = T;
        end
        
        % run predictions, for both birth RFS and existing targets
        function predict(this)
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
            curr_gaussians = repmat(giwComp(0,0,0,0,0,0), 1, n_meas*n_pred);
            
            %update components
            for i = 1:length(this.giw_comps)
                this.giw_comps(i).K = this.giw_comps(i).P*this.H';
                this.giw_comps(i).S = this.H*this.giw_comps(i).K;
                this.giw_comps(i).z = kron(this.H,eye(this.d))*this.giw_comps(i).mu;
                this.giw_comps(i).weight = (1-(1-exp(-this.p_gamma))*this.pd)*this.giw_comps(i).weight;
            end
            
            %update
            for i = 1:n_meas
                for j = 1:n_pred
                    n_points = meas(i).n;
                    S = this.giw_comps(j).S + 1/n_points;
                    inv_S = inv(S);
                    K = this.giw_comps(j).K * inv_S;
                    epsilon = meas(i).center - this.giw_comps(j).z;
                    N = inv_S*epsilon*epsilon';
                    mu = this.giw_comps(j).mu + kron(K,eye(this.d))*epsilon;
                    P = this.giw_comps(j).P - K*S*K';
                    v = this.giw_comps(j).v + n_points;
                    V = this.giw_comps(j).V + N + meas(i).scatter;
                    w = ((exp(-this.p_gamma)*(this.p_gamma)^(n_points)*this.pd)/((this.p_beta^n_points)*((pi^n_points)*n_points*S)^(this.d/2)))...
                        * ((det(this.giw_comps(j).V)^(this.giw_comps(j).v/2))/(det(V)^(v/2)))...
                        * gamma_2d(v/2)/gamma_2d(this.giw_comps(j).v/2)...
                        * this.giw_comps(j).weight;
                end
            end
            
            this.giw_comps = curr_gaussians;
            this.weight_sort_gaussians;
            this.prune;
            this.merge;
        end
        
        function weight_sort_gaussians(this)
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
            
            
            if length(this.giw_comps) > this.max_gaussians             
                this.giw_comps = this.giw_comps(1:this.max_gaussians);
            end
            
            this.recalculate_weights(weightsum_before_prune)
        end
        
        %merging closeby targets
        function merge(this)
            i = 1;
            while i < length(this.giw_comps)
                weightsum = this.giw_comps(i).weight;
                mu_sum = this.giw_comps(i).weight.*this.giw_comps(i).mu;       
                ind = [];
                %check all subsequent gaussians if they are closeby
                for j=(i+1):length(this.giw_comps)
                    proximity = (this.giw_comps(i).mu-this.giw_comps(j).mu)' ...
                        *inv(this.giw_comps(i).P)...
                        *(this.giw_comps(i).mu-this.giw_comps(j).mu);
                    
                    if proximity < this.min_merge_dist
                        weightsum = weightsum + this.giw_comps(j).weight;
                        mu_sum = mu_sum + this.giw_comps(j).weight ...
                            *this.giw_comps(j).mu;
                        ind = [ind j];
                    end                        
                end
                %update mu, P and weight of gaussians(i) if there were
                %closeby gaussians that can be merged
                if ~isempty(ind)
                    mu_sum = mu_sum./weightsum;
                    P_sum = this.giw_comps(i).weight.*this.giw_comps(i).P; 
                    for j=ind
                        P_sum = P_sum + this.giw_comps(j).weight*( ...
                            this.giw_comps(j).P + ( ...
                            mu_sum - this.giw_comps(j).mu) * (...
                            mu_sum - this.giw_comps(j).mu)');
                    end       
                    P_sum = P_sum./weightsum;
                    this.giw_comps(i).mu = mu_sum;
                    this.giw_comps(i).P = P_sum;
                    this.giw_comps(i).weight = weightsum;
                    this.giw_comps(ind) = [];                 
                end
                
                %check if there are still duplicate indices left and
                %reassign their index if need be
                for j=(i+1):length(this.giw_comps)
                    curr_ind = this.giw_comps(i).index;
                    if this.giw_comps(j).index == curr_ind
                        this.index = this.index+1;
                        this.giw_comps(j).index = this.index;
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
            end
        end
        
        function number_of_targets = get_number_of_targets(this)
            number_of_targets = this.number_of_targets;
        end
        
        function best_estimates = get_best_estimates(this)
            n = round(this.number_of_targets);
            best_estimates = [];
            if n > 0
                if n > length(this.giw_comps)
                    n = length(this.giw_comps)
                end
                best_estimates = this.giw_comps(1:n);
            end
        end
        
        function gaussians = get_all_gaussians(this)
            gaussians = this.giw_comps;
        end
    end
end