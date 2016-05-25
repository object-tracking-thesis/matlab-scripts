% Just like MGPgenerator3, but has preprocessing step for each cluster, and
% removes MGPs that are not close enough to measurements. Will also use a
% changed car model.
%
%
%
% Public Methods:
%                              this = MGPgenerator4(N, covMGPgate)
% [gatedMgpHandles, gatedAssignedZ] = generate(clusterZ, predictedState)            
%
% Private Methods:
% [filtClust, massVec1, massVec2, cp] = this.preFilter(clusterZ, lb, ub, gateCov);             
%                  [wVector, lVector] = this.getVectors(predictedState, massVec1, massVec2);            
%
%                              corner = this.getCorner2(predictedState, wVector, lVector);        
%                  [wViewed, lViewed] = this.getViewLen(filtClust, cp, wVector, lVector);             
%
%                           assignedZ = this.selectMeas(filtClust, wVector, lVector, wViewed, lViewed);             
%                          mgpHandles = this.makeMGPs(corner, predictedState, wViewed, lViewed);             
%
%   [gatedMgpHandles, gatedAssignedZ] = this.gateMGPs(assignedZ, mgpHandles);
%

classdef MGPgenerator4 < handle
    properties (Access = public)
        mgpFunCorner1_w; % These will all be generic function handles, dependant on st and K,h,p
        mgpFunCorner1_l;
        
        mgpFunCorner2_w;
        mgpFunCorner2_l;
        
        mgpFunCorner3_w;
        mgpFunCorner3_l;  
        
        mgpFunCorner4_w;
        mgpFunCorner4_l;
         
        mgpNum;
        covMGPgate; % Measurement covariance for each MGP
    end
    
    methods        
        function this = MGPgenerator4(N, covMGPgate)
            % Constructor. Initates the MGPgenerator, where N specifies how
            % many MGPs should be spread out on each side, in addition to
            % the three minimum. 
            this.mgpNum = N;
            this.covMGPgate = covMGPgate;
            % Corner i, width and length MGPs definitions           
            
            this.mgpFunCorner1_w = @(K,h,p) (@(st)... 
                                   [st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + 1/K*h*p*st(6)*cos(st(4) + pi/2);...
                                    st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + 1/K*h*p*st(6)*sin(st(4) + pi/2)]);
                                 
            this.mgpFunCorner1_l = @(K,h,p) (@(st)... 
                                    [st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + 1/K*h*p*st(7)*cos(st(4));...
                                     st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + 1/K*h*p*st(7)*sin(st(4))]);
                                 
            % Corner ii, width and length MGPs definitions

            this.mgpFunCorner2_w = @(K,h,p) (@(st)... 
                                    [st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + st(6)*cos(st(4) + pi/2) + 1/K*h*p*st(6)*cos(st(4) +pi/2 +pi);...
                                     st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + st(6)*sin(st(4) + pi/2) + 1/K*h*p*st(6)*sin(st(4) +pi/2 +pi)]);
                                
            this.mgpFunCorner2_l = @(K,h,p) (@(st)... 
                                    [st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + st(6)*cos(st(4) + pi/2) + 1/K*h*p*st(7)*cos(st(4));...
                                     st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + st(6)*sin(st(4) + pi/2) + 1/K*h*p*st(7)*sin(st(4))]);
                                
            % Corner iii, width and length MGPs definitions

            this.mgpFunCorner3_w = @(K,h,p) (@(st)... 
                                    [st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + st(6)*cos(st(4) + pi/2) + st(7)*cos(st(4)) + 1/K*h*p*st(6)*cos(st(4) + pi/2 + pi);...
                                     st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + st(6)*sin(st(4) + pi/2) + st(7)*sin(st(4)) + 1/K*h*p*st(6)*sin(st(4) + pi/2 + pi)]);
                                
            this.mgpFunCorner3_l = @(K,h,p) (@(st)... 
                                    [st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + st(6)*cos(st(4) + pi/2) + st(7)*cos(st(4)) + 1/K*h*p*st(7)*cos(st(4) + pi);...
                                     st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + st(6)*sin(st(4) + pi/2) + st(7)*sin(st(4)) + 1/K*h*p*st(7)*sin(st(4) + pi)]);
                                
                                
            % Corner iv, width and length MGPs definitions

            this.mgpFunCorner4_w = @(K,h,p) (@(st)... 
                                    [st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + st(7)*cos(st(4)) + 1/K*h*p*st(6)*cos(st(4) + pi/2);...
                                     st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + st(7)*sin(st(4)) + 1/K*h*p*st(6)*sin(st(4) + pi/2)]);

            this.mgpFunCorner4_l = @(K,h,p) (@(st)...
                                    [st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + st(7)*cos(st(4)) + 1/K*h*p*st(7)*cos(st(4) + pi);...
                                     st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + st(7)*sin(st(4)) + 1/K*h*p*st(7)*sin(st(4) + pi)]);

        end
        
        function [gatedMgpHandles, gatedAssignedZ] = generate(this, clusterZ, predictedState)
            % Calls all private functions and returns assigned measurements
            % and their corresponding MPG function handles, for later use
            % in UKF.
            
            % These should not be hardcoded, but submitted by the user
            lb = 0.2;
            ub = 0.5;
            gateCov = 0.2^2;
            
            [filtClust, massVec1, massVec2, cp] = this.preFilter(clusterZ, lb, ub, gateCov);
            
                             [wVector, lVector] = this.getVectors(predictedState, massVec1, massVec2);
            
                                         corner = this.getCorner2(predictedState, wVector, lVector);
            
                             [wViewed, lViewed] = this.getViewLen(filtClust, cp, wVector, lVector);
            
                                      assignedZ = this.selectMeas(filtClust, cp, wVector, lVector, wViewed, lViewed);
            
                                     mgpHandles = this.makeMGPs(corner, predictedState, wViewed, lViewed);
            
              [gatedMgpHandles, gatedAssignedZ] = this.gateMGPs(assignedZ, mgpHandles);
        end
    end
    
    methods (Access = public)
        
        function [filtClust, massVec1, massVec2, cp] = preFilter(~, clusterZ, lb, ub, gateCov)
            
            [~, ~, uOp, ~] = cornerPoint(clusterZ, lb, ub); % 0.2 0.5
            
            c1 = uOp(1); c2 = uOp(2);
            n1 = uOp(3); n2 = uOp(4);
            
            xc = (-n1*c1 + n2*c2);
            yc = (-n2*c1 -n1*c2);
            
            % Formulate all 4 possible L-vectors, to find the ones that are along the data            
            M = mean(clusterZ(:,1:2));
            
            v1 = [1, -n1/n2]; v1 = v1./norm(v1,2);
            v2 = [1,  n2/n1]; v2 = v2./norm(v2,2);
            v3 = -1.*v1;
            v4 = -1.*v2;
            
            p1 = [xc yc] + v1;
            p2 = [xc yc] + v2;
            p3 = [xc yc] + v3;
            p4 = [xc yc] + v4;
            
            d1 = (p1(1) - M(1))^2 + (p1(2) - M(2))^2;
            d2 = (p2(1) - M(1))^2 + (p2(2) - M(2))^2;
            d3 = (p3(1) - M(1))^2 + (p3(2) - M(2))^2;
            d4 = (p4(1) - M(1))^2 + (p4(2) - M(2))^2;
            
            % The massVec* are vectors which are in the direction of the point cloud points, from [xc yc]
            
            if d1 < d3
                massVec1 = v1;
            else
                massVec1 = v3;
            end
            
            if d2 < d4
                massVec2 = v2;
            else
                massVec2 = v4;
            end
            
            %gateCov = 0.2^2;
            filtClust = gateRectangular(clusterZ, [xc, yc], massVec1, massVec2, gateCov);            
            
            cp = [xc,yc];
        end
               
        function [wVector, lVector] = getVectors(~, predictedState, massVec1, massVec2)
            
            vHeading = [cos(predictedState(4)), sin(predictedState(4))];
            % Find which vector is along length & width
            if abs(dot(massVec1, vHeading)) < abs(dot(massVec2,vHeading))
                lVector =  massVec2;
                wVector =  massVec1;
            else
                lVector =  massVec1;
                wVector =  massVec2;
            end
        end
        
        function corner = getCorner2(~, predictedState, wVector, lVector)
            compV = lVector + wVector;
            
            phi = predictedState(4);
    
            R = [cos(phi) -sin(phi); sin(phi) cos(phi)];
            
            % Negate heading            
            compV = compV*R;
            
            angle = mod(atan2(compV(2), compV(1)),2*pi);            
            
            if (0 < angle) && (angle < pi/2) % Quadrant 1
                corner = 1;
            elseif (pi/2 < angle) &&  (angle < pi) % Quadrant 2
                corner = 4;
            elseif (pi < angle) && (angle < 3*pi/2) % Quadrant 3
                corner = 3;
            elseif (3*pi/2 < angle) && (angle < 2*pi) % Quadrant 4
                corner = 2;
            end
        end
        
        function [wViewed, lViewed] = getViewLen(~, filtClust, cp, wVector, lVector)
            % Returns the length of the viewed length and width of the car 
            
            clusterZ = filtClust(:,1:2);

            xc = cp(1);
            yc = cp(2);
                        
            % Project point cloud points onto each line
            top = dot(clusterZ, repmat(wVector, length(clusterZ),1),2);
            widthPoints = [top top] .* repmat(wVector, length(clusterZ), 1);
            
            top = dot(clusterZ, repmat(lVector, length(clusterZ),1),2);
            lengthPoints = [top top] .* repmat(lVector, length(clusterZ), 1);
            
            % Project Corners
            widthCorner  = dot([xc, yc], wVector).*wVector;
            lengthCorner = dot([xc, yc], lVector).*lVector;
                        
            % Calculate euclidian distances between corner and projected points for
            % each projected line
                        
            widthSquaredDist = (repmat(widthCorner(:,1), length(clusterZ), 1) - widthPoints(:,1)).^2 + (repmat(widthCorner(:,2), length(clusterZ), 1) - widthPoints(:,2)).^2;
            
            [~, widthIndexMax] = max(widthSquaredDist);
            
            lengthSquaredDist = (repmat(lengthCorner(:,1), length(clusterZ), 1) - lengthPoints(:,1)).^2 + (repmat(lengthCorner(:,2), length(clusterZ), 1) - lengthPoints(:,2)).^2;
            
            [~, lengthIndexMax] = max(lengthSquaredDist);                      
            
            wViewed = norm(widthCorner - widthPoints(widthIndexMax,:),2);
            
            lViewed = norm(lengthCorner - lengthPoints(lengthIndexMax,:),2);
            
        end
                        
        function assignedZ = selectMeas(this, filtClust, cp, wVector, lVector, wViewed, lViewed)
                                   
            % N = 1; % Test assumption
            N = this.mgpNum;
            
            xc = cp(1);
            yc = cp(2);
            
            storageCP = [xc yc];
            
            storageW = [];
            storageL = [];
                        
            if wViewed > 0.5
                storageW = zeros(N+1,2); % Storage for the practical Width MGPs
                
                for j = 1:N+1
                    storageW(j,:) = [xc, yc] + wViewed/(N+1)*j*wVector;
                end
            end
            
            if lViewed > 0.5
                storageL = zeros(N+1,2); % Storage for the practical Length MGPs                
                
                for j = 1:N+1
                    storageL(j,:) = [xc, yc] + lViewed/(N+1)*j*lVector;
                end
            end
            
            allStorage = [storageL;
                          storageCP;
                          storageW];
            
            % This is where 3sigma gating should occur 
            [idx, dist] = knnsearch(filtClust(:,1:2), allStorage);
            
            Nokeep = dist > 3*sqrt(this.covMGPgate);            
            idxNoKeep = Nokeep.*idx;
            idxNoKeep = idxNoKeep(idxNoKeep > 0);            
            
            % Measurements outside of 3sigma should not be used. Set to NaN
            % so that they can later on be removed together with their
            % associated MGPs.
            filtClust(idxNoKeep, 1:2) = NaN; 
            
            chosenMeasurements = filtClust(idx,1:2);
            
            assignedZ = chosenMeasurements;
            
        end
        
        function [mgpHandles] = makeMGPs(this, corner, predictedState, wViewed, lViewed)
            nrStore = 1+this.mgpNum;
            
            mgpCornerH = cell(1,1);
            
            mgpWidthH = cell(nrStore,1);
            mgpLengthH = cell(nrStore,1);
            
            N = this.mgpNum;
            K = N+1;
            
            if wViewed > 0.5
                w_handle = this.getGeneralHandle(corner, 'w'); % This should return the general function handle 
                % Generate handle for corner 
                mgpCornerH{1} = w_handle(K, 0, 0);
            else
                w_handle = @(K,j,p) [];     % These should be handles to functions that return nothing 
            end
            
            if lViewed > 0.5
                l_handle = this.getGeneralHandle(corner, 'l'); % This should return the general function handle 
                % Generate handle for corner 
                mgpCornerH{1} = l_handle(K, 0, 0);                
            else
                l_handle = @(K,j,p) [];    % These should be handles to functions that return nothing 
            end
                        
            h = 1:K;
            
            for j = h;
                p = wViewed/predictedState(6);
                if p > 1
                    p = 1;
                end
                mgpWidthH{j}  = w_handle(K, j, p);                
                
                p = lViewed/predictedState(7);
                if p > 1
                    p = 1;
                end
                
                mgpLengthH{j}  = l_handle(K, j, p);

            end
            
            % Remove potential empty cells
            mgpLengthH = mgpLengthH(~cellfun(@isempty, mgpLengthH));
            mgpWidthH  = mgpWidthH(~cellfun(@isempty, mgpWidthH));
            
            
            mgpHandles = [mgpLengthH;
                          mgpCornerH; 
                          mgpWidthH];            
                        
        end
        
        function [gatedMgpHandles, gatedAssignedZ] = gateMGPs(~, assignedZ, mgpHandles)
            % Checks if any of the assignedZs have 'NaN' values, indicating
            % that the corresponding MGP doesn't have a measurement within
            % 3sigma range. 
                        
            notNanIndex = ~isnan(assignedZ(:,1));
            
            gatedAssignedZ = assignedZ(notNanIndex, :);
            gatedMgpHandles = mgpHandles(notNanIndex, :);            
            
        end
                        
        %% ====== INTERNAL (hidden) FUNCTIONS ======                
        function [funHandle] = getGeneralHandle(this, corner, side)
            if strcmp(side, 'w')
                switch corner
                    case 1
                        funHandle = this.mgpFunCorner1_w;
                    case 2
                        funHandle = this.mgpFunCorner2_w;
                    case 3
                        funHandle = this.mgpFunCorner3_w;
                    case 4
                        funHandle = this.mgpFunCorner4_w;
                end
                
            elseif strcmp(side, 'l')
                switch corner
                    case 1
                        funHandle = this.mgpFunCorner1_l;
                    case 2
                        funHandle = this.mgpFunCorner2_l;
                    case 3
                        funHandle = this.mgpFunCorner3_l;
                    case 4
                        funHandle = this.mgpFunCorner4_l;
                end
            end
        end
 
    end
end





