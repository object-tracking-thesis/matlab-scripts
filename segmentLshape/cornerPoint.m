function [m1, m2, uOp, varargout] = cornerPoint(pointCloud, varargin)

if length(varargin) == 2
    lb = varargin{1};
    ub = varargin{2};    
else
    lb = 0.15;
    ub = 0.4;
end

Ntg = sortrows(pointCloud,3);
% Remove tailing points

N1 = round(length(Ntg)*lb); % 0.15
N2 = round(length(Ntg)*ub); % 0.4

sortNtg = Ntg(N1:end-N2,1:2);

if nargout == 4    
    varargout(1) = {sortNtg};
end


%% Find clockwise ordering for points

% Get the geometric center point of the pointCloud
m = mean(sortNtg);
% Define perpendicular line
perpLine = [1 -m(1)/m(2)]./norm([1 -m(1)/m(2)],2);

% project points onto perp line, using the formula proj_s-on-v =
% s*dot(v,s) since s is unitvector 

[r,c ] = size(sortNtg);

scalingFactor = sum(sortNtg .* repmat(perpLine,r,1),2); 
projPoints = repmat(scalingFactor,1,c) .* repmat(perpLine,r,1);
%
% Get angle for PC1 & PC2
x = [0 1];
phi = acos(sum(perpLine.*x)/(norm(x,2)*norm(perpLine,2)));
% Rotate PC1 points by phi
R = [cos(phi) -sin(phi); sin(phi) cos(phi)];

projPointsRotated = projPoints*R';

% Add Index before sorting them

projPointsRotated = [projPointsRotated, [1:length(projPointsRotated)]'];

sortedProjPointsRotated = sortrows(projPointsRotated,2);

% idxV represents the mapping from sorted to unsorted points
idxV = sortedProjPointsRotated(:,3);

orderedNtg = sortNtg(idxV,:);

%% Calculate L-shape
uOp = ISED(orderedNtg);

c1 = uOp(1);
c2 = uOp(2);
n1 = uOp(3);
n2 = uOp(4);

xc = (-n1*c1 + n2*c2);
yc = (-n2*c1 -n1*c2);

% Define each vector & normalize it 

V1 = [1 -n1/n2];
V2 = [1 n2/n1];
V1 = V1/norm(V1,2);
V2 = V2/norm(V2,2);

% project points onto perp line, using the formula proj_s-on-v =
% s*dot(v,s)/dot(s,s)
V = {V1 V2};
projPoints = cell(1,2);
projC = cell(1,2);

% Project points 
[r,c ] = size(orderedNtg);

for k = 1:length(V)
    top = sum(orderedNtg .* repmat(V{k},r,1),2);    

    scalingFactor = top;
    projPoints{k} = repmat(scalingFactor,1,c) .* repmat(V{k},r,1);
    % Do the same for (xc,yc)
    top = sum([xc yc] .* V{k});
    scalingFactor = top;
    projC{k} = scalingFactor.* V{k};
end

% Look at length on each projected line 

phi1 = acos(dot(V1, [0 1])); % Get rotation angle 
phi2 = acos(dot(V2, [0 1]));

R1 = [cos(phi1) -sin(phi1); sin(phi1) cos(phi1)];
R2 = [cos(phi2) -sin(phi2); sin(phi2) cos(phi2)];

projPointsR1 = projPoints{1}*R1'; % Apply rotation 
projPointsR2 = projPoints{2}*R2';

projPointsR1 = projPointsR1(abs(projPointsR1) < 1e-5 == 0); % Remove rounding errors 
projPointsR2 = projPointsR2(abs(projPointsR2) < 1e-5 == 0);

projPointsR1 = projPointsR1 - min(projPointsR1); % Set measure from zero 
projPointsR2 = projPointsR2 - min(projPointsR2);

m1 = max(projPointsR1);
m2 = max(projPointsR2);

end