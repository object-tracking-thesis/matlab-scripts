classdef OptimalHypos < handle 
    % Static class that is responsible for calculating optimal hypotheses.
    % Has only one pulic function, generateHypos, which generates
    % hyoptheses.
    
    methods (Access = public)
        
        function optHypos = generateHypos(this, scan, gatingMat, tracks, N)
            % Generates optimal Hypotheses using murty.m and acution.m
            %
            % optHypos = generateHypos(scan, gatingMat, tracks, N)
            %   scan      - scan object containing measurements
            %   gatingMat - gatingMatrix for the measurements and tracks
            %   tracks    - tracks object containing posteriors
            %   N         - Maximal nr of hypotheses wanted to be generated
            %
            betaFA = Model.rho;
            betaNT = Model.spwn;
            x = -1e10; % -inf approximation
            
            matrixFA = diag(betaFA.*ones(size(scan.measId)));
            matrixFA(matrixFA == 0) = x;
            matrixNT = diag(betaNT.*ones(size(scan.measId)));
            matrixNT(matrixNT == 0) = x;
            % gatingMat
            a = this.getTargetAssignmentMat(gatingMat, scan, tracks);
            % matrixFA
            % matrixNT
            % a
            
            assignmentMat = [a, matrixNT, matrixFA];
            
            [tg2m, ~, ~, ~] = Murty(assignmentMat, N);
            
            [~, limitFA] = size([a, matrixNT]);            
            tg2m = tg2m';
            tg2m(tg2m > limitFA) = 0;
            optHypos = tg2m;
        end
        
        function tgAssignmentMat =  getTargetAssignmentMat(this, gatingMat, scan, tracks)
            if isempty(tracks.trackId) || isempty(gatingMat)
                tgAssignmentMat = [];
            else
                [meas, tg] = size(gatingMat);
                tgAssignmentMat = zeros(meas, tg);
                
                for m = 1:meas
                    for t = 1:tg
                        if gatingMat(m, t) == 1 % If this measurement m has been gated with target t
                            tgAssignmentMat(m,t) = this.calcLikelihood(scan.measurements(:,m), tracks.track(t)); % Functioncall to calcLikelihood instead
                        end
                    end
                end
            end
        end
        
        function likelihood = calcLikelihood(~, measurement, track)
            predictedMeasurement = Model.H*track.expectedValue;
            predictedVariance = Model.H*track.covariance*Model.H' + Model.R;
            
            likelihood = mvnpdf(measurement, predictedMeasurement, predictedVariance);
        end
        
    end
end