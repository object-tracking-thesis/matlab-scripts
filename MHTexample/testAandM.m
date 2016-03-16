% Test acution and murty 

function testAandM()

% assMat = [0.0001 0.01 0.1 5;...
%           0.0001 0.01 0.1 5;...
%           0.01   0.1  0.1 5];
% N = 10;
% [tg2m, ~, nrOfSol, ~] = Murty(assMat,N);
% tg2m';

% Test out assignment matrix generation 

betaFA = 0.08;
betaNT = 0.07; 
scan = Scan([[1;2] [3;4]]); % Scan with two measurements 

matrixFA = diag(betaFA*ones(size(scan.measId))); % Create assigment matrix for FAs;
matrixNT = diag(betaNT*ones(size(scan.measId))); % Create assigment matrix för NTs;

% use gatingmatrix to create assigment matrix for existing targets; 

% Let's say we have three targets 
gatingMatrix = [1 1;... % Measurement 1 is gated for T1 and T2 
                1 1];   % Measurement 2 is gated for T2                
                            
% Create targets 
tracks = Tracks;
p1 = Posterior([1 1 1 1]', eye(4));
p2 = Posterior([2.5 2.5 3.5 3.5]', eye(4));
tracks.addTrack(p1);
%tracks.addTrack(p2);

a = getTargetAssigmentMat(gatingMatrix, scan, tracks);
N = 20;
assignmentMat = [a, matrixNT, matrixFA];
[tg2m, ~, ~, ~] = Murty(assignmentMat, N);
disp(assignmentMat)
disp(tg2m')

[~, limitFA] = size([a, matrixNT]);
tg2m = tg2m';
tg2m(tg2m > limitFA) = 0;
disp(tg2m);



end
% pseudo function call 
function tgAssignmentMat =  getTargetAssigmentMat(gatingMat, scan, tracks)
if isempty(tracks.trackId || isempty(gatingMat))
   tgAssignmentMat = [];
 
else
[meas, tg] = size(gatingMat); 
tgAssignmentMat = zeros(size(meas, tg));

for m = 1:meas
    for t = 1:tg
       if gatingMat(m, t) == 1 % If this measurement m has been gated with target t 
           tgAssignmentMat(m,t) = calcLikelihood(scan.measurements(:,m), tracks.track(t)); % Functioncall to calcLikelihood instead
       end        
    end    
end

end
end

function likelihood = calcLikelihood(measurement, track)
    predictedMeasurement = Model.H*track.expectedValue;
    predictedVariance = Model.H*track.covariance*Model.H' + Model.R;
    
    likelihood = mvnpdf(measurement, predictedMeasurement, predictedVariance);
end













