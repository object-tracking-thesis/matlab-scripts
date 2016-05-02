
classdef MGPgenerator < handle
    properties%(Access = private)
        mgpFunX1;
        mgpFunY1;
        
        mgpFunX2;
        mgpFunY2;
        
        mgpFunX3;
        mgpFunY3;
        
        mgpFunX4;
        mgpFunY4;
        
        mgpJac1;
        mgpJac2;
        mgpJac3;
        mgpJac4;
        
        mgpNum;
        
        side1;
        side2;
        side3;
        side4;
    end
    
    methods
        
        function this = MGPgenerator(N)
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
            % ---> 
            % HOWEVER! Duplicates can occur, as such a check will be
            % necessary down the line to only keep unique MGPs.
            syms K h 
            
            % Side 1
            this.mgpFunX1 = x - 0.5*sqrt(W^2 + L^2)*cos(atan(W/L) + phi) + (W*h/K)*cos(phi + pi/2);
            this.mgpFunY1 = y - 0.5*sqrt(W^2 + L^2)*sin(atan(W/L) + phi) + (W*h/K)*sin(phi + pi/2);
            
            % Side 2
            this.mgpFunX2 = x - 0.5*sqrt(W^2 + L^2)*cos(atan(W/L) + phi) + (W)*cos(phi + pi/2) + (L*h/K)*cos(phi);
            this.mgpFunY2 = y - 0.5*sqrt(W^2 + L^2)*sin(atan(W/L) + phi) + (W)*sin(phi + pi/2) + (L*h/K)*sin(phi);
            
            % Side 3
            this.mgpFunX3 = x - 0.5*sqrt(W^2 + L^2)*cos(atan(W/L) + phi) + (W*h/K)*cos(phi + pi/2) + (L)*cos(phi);
            this.mgpFunY3 = y - 0.5*sqrt(W^2 + L^2)*sin(atan(W/L) + phi) + (W*h/K)*sin(phi + pi/2) + (L)*sin(phi);
            
            % Side 4
            this.mgpFunX4 = x - 0.5*sqrt(W^2 + L^2)*cos(atan(W/L) + phi) + (L*h/K)*cos(phi);
            this.mgpFunY4 = y - 0.5*sqrt(W^2 + L^2)*sin(atan(W/L) + phi) + (L*h/K)*sin(phi);
            
            this.mgpJac1 = this.formJacobian(this.mgpFunX1, this.mgpFunY1);
            this.mgpJac2 = this.formJacobian(this.mgpFunX2, this.mgpFunY2);
            this.mgpJac3 = this.formJacobian(this.mgpFunX3, this.mgpFunY3);
            this.mgpJac4 = this.formJacobian(this.mgpFunX4, this.mgpFunY4);  
            
            % Define MGP functions for each car middle side 

            %    _<-3_
            % ^ |     | ^
            % | |  ?  | |
            % 2 |  x  | 4
            %   |     |
            %   |_____|
            %    <-1
            %
            %       1    2    3    4      5         6    7
            % st = [x_k, y_k, v_k, phi_k, phiDot_k, w_k, l_k]';

            this.side1 = @(st) ...
                [...
                st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + 0.5*st(6)*cos(st(4) + pi/2);...
                st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + 0.5*st(6)*sin(st(4) + pi/2)
                ];

            this.side2 = @(st) ... 
                [...
                st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + st(6)*cos(st(4) + pi/2) + 0.5*st(7)*cos(st(4));...
                st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + st(6)*sin(st(4) + pi/2) + 0.5*st(7)*sin(st(4))
                ];

            this.side3 = @(st) ...
                [...
                st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + 0.5*st(6)*cos(st(4) + pi/2) + st(7)*cos(st(4));...
                st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + 0.5*st(6)*sin(st(4) + pi/2) + st(7)*sin(st(4))
                ];

            this.side4 = @(st)...
                [...
                st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + 0.5*st(7)*cos(st(4));...
                st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + 0.5*st(7)*sin(st(4))
                ];
            
            
        end
        
        function vehicleCase = decideCase(this,clusterZ, predictedState)
            % NOTE!! clusterZ must be xy only, not z
            % Takes the clusterZ, performs an L-shape fitting to find the
            % corner point, and with this we construct a rectangle car around the
            % predicted state, and then chose to keep the two sides that
            % are nearest to the corner point. If either of the length 
            % measurements is below a certain threshold then only one of
            % the sides is chosen. < 0.5 m seems good.            
            
            % Get cornerpoint             
            [m1,m2,uOp] = cornerPoint(clusterZ); 
            
            
            c1 = uOp(1); c2 = uOp(2);
            n1 = uOp(3); n2 = uOp(4);

            xc = (-n1*c1 + n2*c2); % Define cornerpoint 
            yc = (-n2*c1 -n1*c2);
            
            cPoint = [xc, yc];
            
            % calculate middle MGP for each side
            allSides = [[this.side1(predictedState)]';...
                        [this.side2(predictedState)]';...
                        [this.side3(predictedState)]';... 
                        [this.side4(predictedState)]'];            
            % If either measure less than 0.5m, then only one side should
            % get MGPs. To find that side, we take the MGP that is closest
            % to the middle value. Otherwise, we take the two middle MGPs
            % that are closest to the cornerpoint
            if (m1 < 0.5) || (m2 < 0.5)
                pcMean = mean(clusterZ);
                D = zeros(4,1);
                for j = 1:4
                    D(j) = sqrt((allSides(j,1) - pcMean(1))^2 +... 
                                (allSides(j,2) - pcMean(2))^2);
                end
                
                [~, idx] = min(D);
                vehicleCase = idx;
            else
                D = zeros(1,4);
                for j = 1:4
                    D(j) = sqrt((allSides(j,1) - cPoint(1))^2 +... 
                                (allSides(j,2) - cPoint(2))^2);
                end                
                [~, idx1] = min(D);                
                D(idx1) = inf;
                [~, idx2] = min(D);
                vehicleCase = [idx1 idx2];
            end
        end
        
        function [mgps, jacobians] = constructMGPs(this, vehicleCase, predictedState)
            % Using the vehicle case, as well as the predicted state, the
            % function constructs the MGPs and their jacobians based on the
            % parameters in the predicted state. As such, eval calls are
            % made on the symbolic functions defined in the constructor.
            
            K = (this.mgpNum+1);
            totalNrMgps = length(vehicleCase)*(K+1);
            
            mgps = zeros(totalNrMgps,2);
            jacobians = zeros(2,7,totalNrMgps);
            
            index = 1;
            
            for j = vehicleCase
                % Get general symbolic functions for MGPs along the car side j
                [funXsym, funYsym, symJacob] = this.getSymFunAndJacob(j);
                
                for h = 0:K                                      
                    xh = this.evaluateFunction(funXsym, predictedState, K, h);
                    yh = this.evaluateFunction(funYsym, predictedState, K, h);
                    Jh = this.evaluateJacobian(symJacob, predictedState, K, h);
                    
                    mgps(index,:) = [xh, yh];                    
                    jacobians(:,:,index) = Jh;
                    
                    index = index+1;
                end                
            end
            
            % Check and see if there are any duplicates (this can happen due to how the functions are defined)
            [~, idx] = unique(mgps, 'rows'); % Find index of unique rows (i.e. MGPs)
            mgps = mgps(idx,:);
            jacobians = jacobians(:,:,idx);
            
        end
        
        function [assignedMgps, assignedJacobians, assignedZ] = assignMgps(~, mgps, jacobians, clusterZ)
            % Takes each generated MGP, and assigns them the
            % closest measurement (euclidian distance).
            % MGPs can't share points, and as such assignment is done
            % iteratively, where the MGP who is closest to a point gets it,
            % then the next MGP is assigned and so on.
            clusterZ = clusterZ(:,1:2);
            
            assignedMgps = ones(size(mgps));
            assignedJacobians = ones(size(jacobians));
            assignedZ = ones(size(assignedMgps));
            
            nrOfMgps = size(assignedMgps,1);
            
            for j = 1:nrOfMgps % nr of rows, i.e. MGPS
                
                % For each mgp, find the corresponding distance and idx in
                % clusterZ that is shortest
                [idx, distance] = knnsearch(clusterZ, mgps);

                [~,I] = min(distance); % Find index of shortest distance                                                
                                
                assignedMgps(j,:) = mgps(I,:); % Keep shortest distance mgp
                assignedJacobians(:,:,j) = jacobians(:,:,I); % Keep shortast distance mgp jacobian
                assignedZ(j,:) = clusterZ(idx(I), :); % Keep corresponding measurement to shortest distance 
                
                mgps(I,:) = [];
                jacobians(:,:,I) = [];
                clusterZ(idx(I), :) = [];                
                
            end
        end
        
        %% ====== API FUNCTIONS ======
        
        function symbolicJacob = formJacobian(~, funX, funY)
            syms x y v phi phiDot W L
            symbolicJacob = [diff(funX, x),...
                             diff(funX, y),...
                             diff(funX, v),...
                             diff(funX, phi),...
                             diff(funX, phiDot),...
                             diff(funX, W),...
                             diff(funX, L);...
                             % second row
                             diff(funY, x),...
                             diff(funY, y),...
                             diff(funY, v),...
                             diff(funY, phi),...
                             diff(funY, phiDot),...
                             diff(funY, W),...
                             diff(funY, L)];
        end
        
        function numericJacob = evaluateJacobian(~, symbolicJacobian, state, K, h) %#ok<INUSD>
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
        
        function numericFunction = evaluateFunction(~, symbolicFunction, state, K, h) %#ok<INUSD>
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
        
        function [funXsym, funYsym, symJacob] = getSymFunAndJacob(this, vehicleCase)
            switch vehicleCase
                case 1
                    funXsym = this.mgpFunX1;
                    funYsym = this.mgpFunY1;
                    symJacob = this.mgpJac1;
                case 2
                    funXsym = this.mgpFunX2;
                    funYsym = this.mgpFunY2;
                    symJacob = this.mgpJac2;
                case 3
                    funXsym = this.mgpFunX3;
                    funYsym = this.mgpFunY3;
                    symJacob = this.mgpJac3;
                case 4
                    funXsym = this.mgpFunX4;
                    funYsym = this.mgpFunY4;
                    symJacob = this.mgpJac4;
            end
        end
    end
end




























