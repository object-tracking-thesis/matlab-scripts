% Class API functions, and order in which they should be called. 
%
%
%                         this = MGPgenerator2(N)
%                       corner = getCarCorner(clusterZ, predictedState)
%           [wViewed, lViewed] = getViewedLengths(clusterZ, predictedState)
% [orderedMgps, orderedJacobs] = constructMGPs(corner, predictedState, w_viewed, l_viewed)
%                    assignedZ = assignMgps(clusterZ, predictedState)
%
% ============= TODO ============
% The way that the class has changed from v1 has made it so that N cant be
% equal to 0. This should be fixed, since as of now the lowest possible nr
% of MGPS is 3+2*1 = 5.
%
classdef MGPgenerator2 < handle
    properties%(Access = private)
        mgpFunCorner1_w;
        mgpFunCorner1_l;
        
        mgpFunCorner2_w;
        mgpFunCorner2_l;
        
        mgpFunCorner3_w;
        mgpFunCorner3_l;  
        
        mgpFunCorner4_w;
        mgpFunCorner4_l;
        
        
        mgpJac1_w;
        mgpJac1_l;
        
        mgpJac2_w;
        mgpJac2_l;
        
        mgpJac3_w;
        mgpJac3_l;
        
        mgpJac4_w;
        mgpJac4_l;
        
        mgpNum;
    end
    
    methods
        
        function this = MGPgenerator2(N)
            this.mgpNum = N;
            % Needs to evaluate and initialize inital functions and their
            % jacobians, where N sets how many mgp's we want between on 
            % each side. In other words, at least 2 or 3 MGP's are garantued, since
            % we will either measure 1 or 2 sides of the car. 
            % As such, we define symbolic jacobians and functions for
            % each general mgp, and allow for an eval call when state is
            % specified, as well as which mgp we are after
            
            % These are state parameters,
            syms x y v phi phiDot W L
            
            % These are parameters for constructing the spread of the MGPs,
            % i.e. K = N+1 & h = 0:K, this will allow us to generate 
            % uniformly spread out mgps on each of the cars sides.
            % p will stand for percentage of side seen 
            % ---> 
            % HOWEVER! Duplicates can occur, as such a check will be
            % necessary down the line to only keep unique MGPs.            
            syms K h p 
            
            % Corner i, width and length MGPs definitions
            this.mgpFunCorner1_w = [ (x - 0.5*sqrt(W^2 + L^2)*cos(atan(W/L) + phi) + 1/K*h*p*W * cos(phi + pi/2)),...
                                     (y - 0.5*sqrt(W^2 + L^2)*sin(atan(W/L) + phi) + 1/K*h*p*W * sin(phi + pi/2))];
            
            this.mgpFunCorner1_l = [ (x - 0.5*sqrt(W^2 + L^2)*cos(atan(W/L) + phi) + 1/K*h*p*L * cos(phi)),...
                                     (y - 0.5*sqrt(W^2 + L^2)*sin(atan(W/L) + phi) + 1/K*h*p*L * sin(phi))];
            
            % Corner ii, width and length MGPs definitions
            this.mgpFunCorner2_w = [(x - 0.5*sqrt(W^2 + L^2)*cos(atan(W/L) + phi) + W*cos(phi + pi/2) + 1/K*h*p*W * cos(phi + pi/2 + pi)),...
                                    (y - 0.5*sqrt(W^2 + L^2)*sin(atan(W/L) + phi) + W*sin(phi + pi/2) + 1/K*h*p*W * sin(phi + pi/2 + pi))];

            this.mgpFunCorner2_l = [(x - 0.5*sqrt(W^2 + L^2)*cos(atan(W/L) + phi) + W*cos(phi + pi/2) + 1/K*h*p*L * cos(phi)),...
                                    (y - 0.5*sqrt(W^2 + L^2)*sin(atan(W/L) + phi) + W*sin(phi + pi/2) + 1/K*h*p*L * sin(phi))];
                        
            % Corner iii, width and length MGPs definitions
            this.mgpFunCorner3_w = [(x - 0.5*sqrt(W^2 + L^2)*cos(atan(W/L) + phi) + W*cos(phi + pi/2) + L*cos(phi) + 1/K*h*p*W * cos(phi + pi/2 + pi)),...
                                    (y - 0.5*sqrt(W^2 + L^2)*sin(atan(W/L) + phi) + W*sin(phi + pi/2) + L*sin(phi) + 1/K*h*p*W * sin(phi + pi/2 + pi))];

            this.mgpFunCorner3_l = [(x - 0.5*sqrt(W^2 + L^2)*cos(atan(W/L) + phi) + W*cos(phi + pi/2) + L*cos(phi) + 1/K*h*p*L * cos(phi + pi)),...
                                    (y - 0.5*sqrt(W^2 + L^2)*sin(atan(W/L) + phi) + W*sin(phi + pi/2) + L*sin(phi) + 1/K*h*p*L * sin(phi + pi))];

            % Corner iv, width and length MGPs definitions
            this.mgpFunCorner4_w = [(x - 0.5*sqrt(W^2 + L^2)*cos(atan(W/L) + phi) + L*cos(phi) + 1/K*h*p*W * cos(phi + pi/2)),...
                                    (y - 0.5*sqrt(W^2 + L^2)*sin(atan(W/L) + phi) + L*sin(phi) + 1/K*h*p*W * sin(phi + pi/2))];

            this.mgpFunCorner4_l = [(x - 0.5*sqrt(W^2 + L^2)*cos(atan(W/L) + phi) + L*cos(phi) + 1/K*h*p*L * cos(phi + pi)),...
                                    (y - 0.5*sqrt(W^2 + L^2)*sin(atan(W/L) + phi) + L*sin(phi) + 1/K*h*p*L * sin(phi + pi))];
            
            % Get symbolic jacobians for each case
            this.mgpJac1_w = jacobian(this.mgpFunCorner1_w, [x y v phi phiDot W L]);
            this.mgpJac1_l = jacobian(this.mgpFunCorner1_l, [x y v phi phiDot W L]);
            
            this.mgpJac2_w = jacobian(this.mgpFunCorner2_w, [x y v phi phiDot W L]);
            this.mgpJac2_l = jacobian(this.mgpFunCorner2_l, [x y v phi phiDot W L]);
            
            this.mgpJac3_w = jacobian(this.mgpFunCorner3_w, [x y v phi phiDot W L]);
            this.mgpJac3_l = jacobian(this.mgpFunCorner3_l, [x y v phi phiDot W L]);
            
            this.mgpJac4_w = jacobian(this.mgpFunCorner4_w, [x y v phi phiDot W L]);
            this.mgpJac4_l = jacobian(this.mgpFunCorner4_l, [x y v phi phiDot W L]);
                
        end
        
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
            p4 = [0 1];
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
        
        function [orderedMgps, orderedJacobs] = constructMGPs(this, corner, predictedState, w_viewed, l_viewed)
            nrStore = 1+this.mgpNum;
            
            mgps_w = nan*zeros(nrStore,2);
            mgps_l = nan*zeros(nrStore,2);
            
            jac_w = zeros(2,7,nrStore)*nan; % Just so that we can keep track of shit
            jac_l = zeros(2,7,nrStore)*nan;
            
            N = this.mgpNum;
            K = N+1;
            
            if w_viewed > 0.5
                [symFun_w, symJac_w] = this.getSymFunAndJacob(corner, 'w');
                mgp_corner = this.evaluateFunction(symFun_w, predictedState, K, 0, 0);
                jac_corner = this.evaluateJacobian(symJac_w, predictedState, K, 0, 0);
            else
                symFun_w = @(x) [];
                symJac_w = @(x) [];
            end
            
            if l_viewed > 0.5
                [symFun_l, symJac_l] = this.getSymFunAndJacob(corner, 'l');
                mgp_corner = this.evaluateFunction(symFun_l, predictedState, K, 0, 0);
                jac_corner = this.evaluateJacobian(symJac_l, predictedState, K, 0, 0);
            else
                symFun_l = @(x) [];
                symJac_l = @(x) [];
            end
                        
            h = 1:K;
            
            for j = h;
                p = w_viewed/predictedState(6);
                if p > 1
                    p = 1;
                end
                mgps_w(j,:)  = this.evaluateFunction(symFun_w, predictedState, K, j, p);                
                jac_w(:,:,j) = this.evaluateJacobian(symJac_w, predictedState, K, j, p);
                
                p = l_viewed/predictedState(7);
                if p > 1
                    p = 1;
                end
                mgps_l(j,:)  = this.evaluateFunction(symFun_l, predictedState, K, j, p);
                jac_l(:,:,j) = this.evaluateJacobian(symJac_l, predictedState, K, j, p);                
            end

            orderedMgps = [mgps_l;
                           mgp_corner; 
                           mgps_w];
            orderedJacobs = cat(3, jac_l, jac_corner, jac_w);
            
            % -----
            % DO WE EVEN NEED THIS ANYMORE?????
            % -----
            % Find unique MGPS
            [~, idx] = uniquetol(orderedMgps, 'ByRows',true); % Find index of unique rows (i.e. MGPs)
            idx = sort(idx);
            orderedMgps = orderedMgps(idx,:);
            orderedJacobs = orderedJacobs(:,:,idx);
            
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
        
        function numericJacob = evaluateJacobian(~, symbolicJacobian, state, K, h, p) %#ok<INUSD>
            % Evaluates symbolic jacobian with the values found in state
            % and K & h
            
            x = state(1);       %#ok<NASGU>
            y = state(2);       %#ok<NASGU>
            v = state(3);       %#ok<NASGU>
            phi = state(4);     %#ok<NASGU>
            phiDot = state(5);  %#ok<NASGU>
            W = state(6);       %#ok<NASGU>
            L = state(7);       %#ok<NASGU>
            
            numericJacob = eval(symbolicJacobian);                        
        end
        
        function numericFunction = evaluateFunction(~, symbolicFunction, state, K, h, p) %#ok<INUSD>
            % Evaluates symbolic jacobian with the values found in state
            % and K & h
            
            x = state(1);       %#ok<NASGU>
            y = state(2);       %#ok<NASGU>
            v = state(3);       %#ok<NASGU>
            phi = state(4);     %#ok<NASGU>
            phiDot = state(5);  %#ok<NASGU>
            W = state(6);       %#ok<NASGU>
            L = state(7);       %#ok<NASGU>
            
            numericFunction = eval(symbolicFunction);                           
        end
        
        function [funSym, symJacob] = getSymFunAndJacob(this, corner, side)
            if strcmp(side, 'w')
                switch corner
                    case 1
                        funSym   = this.mgpFunCorner1_w;
                        symJacob = this.mgpJac1_w;
                    case 2
                        funSym   = this.mgpFunCorner2_w;
                        symJacob = this.mgpJac2_w;
                    case 3
                        funSym   = this.mgpFunCorner3_w;
                        symJacob = this.mgpJac3_w;
                    case 4
                        funSym   = this.mgpFunCorner4_w;
                        symJacob = this.mgpJac4_w;
                end
                
            elseif strcmp(side, 'l')
                switch corner
                    case 1
                        funSym   = this.mgpFunCorner1_l;
                        symJacob = this.mgpJac1_l;
                    case 2
                        funSym   = this.mgpFunCorner2_l;
                        symJacob = this.mgpJac2_l;
                    case 3
                        funSym   = this.mgpFunCorner3_l;
                        symJacob = this.mgpJac3_l;
                    case 4
                        funSym   = this.mgpFunCorner4_l;
                        symJacob = this.mgpJac4_l;
                end
            end
        end
 
    end
end





