classdef RectPHDfilter < handle
    properties
        rectangles = [];
        birth_rectangles = [];
        d = 2;
        ps = 0.98;
        min_survival_weight = 0.0001;
        min_merge_dist = 0.05;
        max_comps = 50;
        number_of_targets = 0;
        index = 1;
        kk = 0
    end
    
    methods(Access = public)
        function this = RectPHDfilter()
            
        end        
        
        % set the birth RFS
        function set_birth_rfs(this, means)
            this.rectangles = [];
            for i = 1:length(means)
                birth_target = RectTarget();
                birth_target.init(means{i}, 0, this.index, 1);
                this.index = this.index+1;
                this.birth_rectangles = [this.birth_rectangles birth_target];
            end 
        end
        
        % run predictions, for existing targets
        function predict(this)
            %new birth_rfs indices
            for i = 1:length(this.birth_rectangles)
                this.birth_rectangles(i).index = this.index;
                this.birth_rectangles(i).predict();
                this.index = this.index+1;
            end
            %prediction for all previous targets
            for i = 1:length(this.rectangles)
                this.rectangles(i).weight = this.ps*this.rectangles(i).weight;
                this.rectangles(i).predict();
            end
            this.rectangles = [this.rectangles this.birth_rectangles];
        end
        
        %update with current measurement Z and isObject vector that
        %contains the NN predictions for each of the measurements with 1
        %being clutter, 2 car, 3 cycle, 4 pedestrian, 5 pedestrian groups
        function update(this, meas)
            %set pd dynamically depending on whether a measurement was
            %received inside the gate of that target and update the weight
            rectangles_no_gating = [];
%             for i = 1:length(this.rectangles)
%                 gating = 0;
%                 for j = 1:length(meas)
%                     gating = this.rectangles(i).gating(meas(j));
%                     if gating
%                         break;
%                     end
%                 end
%                 %when there was no measurement within the gate              
%                 if ~gating
%                     this.rectangles(i).updateWeightNoGating;
%                     ellipse_copy = this.rectangles(i).copy;
%                     rectangles_no_gating = [rectangles_no_gating ellipse_copy];                    
%                 end
%             end
            
            %preallocation for all targets to be updated
            n_meas = length(meas);
            n_pred = length(this.rectangles);
            
            %TODO preallocate here
            new_rectangles = [];
            for i = 1:n_meas
                for j = 1:n_pred
                    ellipse_copy = this.rectangles(j).copy;
                    new_rectangles = [new_rectangles ellipse_copy];
                end
            end            
                   
            %likelihood update
            for i = 1:n_meas      
                weightsum = 0;
                for j = 1:n_pred       
                    this.rectangles(j).weight
                    w = this.rectangles(j).weight * new_rectangles(j+((i-1)*n_pred)).calcLikelihood(meas(i))             
                    new_rectangles(j+((i-1)*n_pred)).weight = w;
                    weightsum = weightsum+w;
                end
                weightsum
                %normalizing the weights with the weightsum for the entire measurement  
                for j=1:n_pred
                    j
                    new_rectangles(j+((i-1)*n_pred)).weight
                    new_rectangles(j+((i-1)*n_pred)).weight = ...
                        new_rectangles(j+((i-1)*n_pred)).weight/(this.kk+weightsum);
                    new_rectangles(j+((i-1)*n_pred)).weight
                end 
            end
            
            this.rectangles = new_rectangles;
            this.weight_sort_components;
            this.prune;
            
            %run the actual update on all new rectangles that actually
            %survived the pruning
            for i = 1:length(this.rectangles)
                this.rectangles(i).update();
            end            
            
            this.rectangles = [this.rectangles  rectangles_no_gating];
            
            %keep only the max_comps best components
            if length(this.rectangles) > this.max_comps             
                this.rectangles = this.rectangles(1:this.max_comps);
            end      
        end
        
        function weight_sort_components(this)
            weightvec = zeros(length(this.rectangles),1);
            for i=1:length(this.rectangles)
                weightvec(i) = this.rectangles(i).weight;
            end
            [~, ind] = sort(weightvec(:),'descend');               
            this.rectangles = this.rectangles(ind);
        end
        
        %prune all gaussians that have a weight below a globally defined
        %threshold
        function prune(this)
            weightvec = zeros(length(this.rectangles),1);
            to_delete = [];
            for i=1:length(this.rectangles)
                weightvec(i) = this.rectangles(i).weight;
                if this.rectangles(i).weight < this.min_survival_weight
                    to_delete = [to_delete i];
                end
            end
            this.rectangles(to_delete) = [];
            weightvec
            weightsum_before_prune = sum(weightvec);
            this.number_of_targets = weightsum_before_prune;                                   
        end
             
        function number_of_targets = get_number_of_targets(this)
            number_of_targets = this.number_of_targets;
        end
        
        function best_estimates = get_best_estimates(this)            
            this.weight_sort_components;
            n = round(this.number_of_targets);
            best_estimates = [];
            if n > 0
                if n > length(this.rectangles)
                    n = length(this.rectangles);
                end
                best_estimates = this.rectangles(1:n);
            end
        end
        
        function gaussians = get_all_gaussians(this)
            gaussians = this.rectangles;
        end
    end
end