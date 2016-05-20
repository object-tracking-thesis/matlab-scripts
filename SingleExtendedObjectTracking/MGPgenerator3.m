% Just like MGPgenerator2, but this version is instead intended to work
% with UKF and not EKF. The difference is that instead of returning values
% of h(state) for each mgp, and their jacobians, this function returns
% function-handles for each mgp-function. 
%
% Public Methods:
%                    this = MGPgenerator3(N)
% [mgpHandles, assignedZ] = generate(clusterZ, predictedState)            
%
% Private Methods:
%             corner = getCarCorner(clusterZ, predictedState)
%       [mgpHandles] = constructMGPs(corner, predictedState, w_viewed, l_viewed)
%          assignedZ = assignMgps(clusterZ, predictedState)
% [wViewed, lViewed] = getViewedLengths(clusterZ, predictedState)
%

classdef MGPgenerator3 < handle
    properties (Access = private)
        mgpFunCorner1_w; % These will all be generic function handles, dependant on st and K,h,p
        mgpFunCorner1_l;
        
        mgpFunCorner2_w;
        mgpFunCorner2_l;
        
        mgpFunCorner3_w;
        mgpFunCorner3_l;  
        
        mgpFunCorner4_w;
        mgpFunCorner4_l;
         
        mgpNum;
    end
    
    methods        
        function this = MGPgenerator3(N)
            % Constructor. Initates the MGPgenerator, where N specifies how
            % many MGPs should be spread out on each side, in addition to
            % the three minimum. 
            this.mgpNum = N;

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
        
        function [mgpHandles, assignedZ] = generate(this, clusterZ, predictedState)
            % Calls all private functions and returns assigned measurements
            % and their corresponding MPG function handles, for later use
            % in UKF.
            corner    = this.getCarCorner(clusterZ, predictedState);
            [wV, lV]  = this.getViewedLengths(clusterZ, predictedState);
            assignedZ = this.assignMgps(clusterZ, predictedState);
            
            [mgpHandles] = this.constructMGPs(corner, predictedState, wV, lV);
        end
    end
    
    methods (Access = private)
        function corner = getCarCorner(~,clusterZ, predictedState) % clusterZ needs Z-axis as well
            % Car corner definition
            % 
            % ii                   iii
            %  *--------------------*
            %  |                    |
            %  |    ----->          |
            %  |                    |
            %  *--------------------*
            %  i                   iv
            
            % Define rectangle for car
            p1 = [0 0];
            p2 = [0 0.4];
            p3 = [1 0.4];
            p4 = [1 0];
            rect = [p1; p2; p3; p4];
            
            
            rot = @(phi) [cos(phi) -sin(phi);
                          sin(phi)  cos(phi)];
            
            trans = @(deltaP) repmat(deltaP, 4,1);
            
            phi = predictedState(4); % Predicted heading
            
            deltaP = mean(clusterZ(:,1:2)); % Mean value of point cloud L-shape
            
            modRect = rect*rot(phi)' + trans(deltaP);
            
            % cC is the modified rectangle, cP is the cornerPoint
            distCorner = @(cC, cP) [(sqrt((cC(1,1) - cP(1))^2 + (cC(1,2) - cP(2))^2));...
                                    (sqrt((cC(2,1) - cP(1))^2 + (cC(2,2) - cP(2))^2));...
                                    (sqrt((cC(3,1) - cP(1))^2 + (cC(3,2) - cP(2))^2));...
                                    (sqrt((cC(4,1) - cP(1))^2 + (cC(4,2) - cP(2))^2))];
            
            [~,~,uOp] = cornerPoint(clusterZ);
            
            c1 = uOp(1); c2 = uOp(2);
            n1 = uOp(3); n2 = uOp(4);
            
            xc = (-n1*c1 + n2*c2);
            yc = (-n2*c1 -n1*c2);
            
            [~, corner] = min(distCorner(modRect, [xc yc]));
            
        end
        
        function [mgpHandles] = constructMGPs(this, corner, predictedState, wViewed, lViewed)
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
        
        function assignedZ = assignMgps(this, clusterZ, predictedState)
            
            [~, ~, uOp] = cornerPoint(clusterZ);
            
            c1 = uOp(1); c2 = uOp(2);
            n1 = uOp(3); n2 = uOp(4);
            
            xc = (-n1*c1 + n2*c2);
            yc = (-n2*c1 -n1*c2);            
            
            M = mean(clusterZ(:,1:2));
            
            % Formulate all 4 possible L-vectors, to find the ones that are along the data
            
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
            
            
            % Define the heading vector
            vHeading = [cos(predictedState(4)), sin(predictedState(4))];
            
            % Find the vectors aligned with length & width of L-shape
            % Orthogonal vectors have dot product ~ 0, i.e. min dot product must be the
            % width vector, since it is orthogonal to the heading vector.
            if abs(dot(massVec1, vHeading)) < abs(dot(massVec2,vHeading))
                vLength =  massVec2;
                vWidth  =  massVec1;
            else
                vLength =  massVec1;
                vWidth  =  massVec2;
            end
            
            % N = 1; % Test assumption
            N = this.mgpNum;
            % Get the viewed lengths of each side
            
            [wViewed, lViewed] = this.getViewedLengths(clusterZ, predictedState);
            
            %wViewed = 1.6016;
            %lViewed = 4.2166;
            
            storageCP = [xc, yc];
            
            storageW = [];
            storageL = [];
            
            
            if wViewed > 0.5
                storageW = zeros(N+1,2); % Storage for the practical Width MGPs
                
                for j = 1:N+1
                    storageW(j,:) = [xc, yc] + wViewed/(N+1)*j*vWidth;
                end
            end
            
            if lViewed > 0.5
                storageL = zeros(N+1,2); % Storage for the practical Length MGPs                
                
                for j = 1:N+1
                    storageL(j,:) = [xc, yc] + lViewed/(N+1)*j*vLength;
                end
            end
            
            allStorage = [storageL;
                storageCP;
                storageW];
            
            IDX = knnsearch(clusterZ(:,1:2), allStorage);
            
            chosenMeasurements = clusterZ(IDX,1:2);
            
            assignedZ = chosenMeasurements;
            
        end
        
        function [wViewed, lViewed] = getViewedLengths(~, clusterZ, predictedState)
            % Returns the length of the viewed length and width of the car 
            
            [~, ~, uOp] = cornerPoint(clusterZ);                        
            clusterZ = clusterZ(:,1:2);
            c1 = uOp(1); c2 = uOp(2);
            n1 = uOp(3); n2 = uOp(4);
            
            xc = (-n1*c1 + n2*c2);
            yc = (-n2*c1 -n1*c2);
            
            % Define the two vectors spanning the L-shape, normalize them
            v1 = [1, -n1/n2]; v1 = v1./norm(v1,2);
            v2 = [1,  n2/n1]; v2 = v2./norm(v2,2);
            
            % Define the heading vector
            vHeading = [cos(predictedState(4)), sin(predictedState(4))];
            
            % Find the vectors aligned with length & width of L-shape
            % Orthogonal vectors have dot product ~ 0, i.e. min dot product must be the
            % width vector, since it is orthogonal to the heading vector.
            if abs(dot(v1, vHeading)) < abs(dot(v2,vHeading))
                vLength =  v2;
                vWidth  =  v1;
            else
                vLength =  v1;
                vWidth  =  v2;
            end
            
            % Project point cloud points onto each line
            top = dot(clusterZ, repmat(vWidth, length(clusterZ),1),2);
            widthPoints = [top top] .* repmat(vWidth, length(clusterZ), 1);
            
            top = dot(clusterZ, repmat(vLength, length(clusterZ),1),2);
            lengthPoints = [top top] .* repmat(vLength, length(clusterZ), 1);
            
            % Project Corners
            widthCorner  = dot([xc, yc], vWidth).*vWidth;
            lengthCorner = dot([xc, yc], vLength).*vLength;
            
            % Calculate euclidian distances between corner and projected points for
            % each projected line
                        
            widthSquaredDist = (repmat(widthCorner(:,1), length(clusterZ), 1) - widthPoints(:,1)).^2 + (repmat(widthCorner(:,2), length(clusterZ), 1) - widthPoints(:,2)).^2;
            
            [~, widthIndexMax] = max(widthSquaredDist);
            
            lengthSquaredDist = (repmat(lengthCorner(:,1), length(clusterZ), 1) - lengthPoints(:,1)).^2 + (repmat(lengthCorner(:,2), length(clusterZ), 1) - lengthPoints(:,2)).^2;
            
            [~, lengthIndexMax] = max(lengthSquaredDist);            
            
            wViewed = norm([xc, yc] - clusterZ(widthIndexMax,:),2);
            
            lViewed = norm([xc, yc] - clusterZ(lengthIndexMax,:),2);
            
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





