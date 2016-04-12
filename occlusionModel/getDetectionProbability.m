function Pd = getDetectionProbability(pointClouds)
% Given a struct containing pointClouds, the function returns the detection
% probability Pd for each cloud.

% Get the number of clouds present in pointClouds
nrClouds = length(pointClouds);

% Allocate space for the two unitvectors that span each cone of each
% pointCloud.
uv1 = zeros(nrClouds,2);
uv2 = zeros(nrClouds,2);

% Allocate space for the targetVector, which is the vector to the geometric
% mean of each pointCloud.
tgVector = zeros(nrClouds,2);
% Allocate space for the corresponing distance from origin to the target
% vector.
distance = zeros(nrClouds,1);
% Allocate space for storing the angle between the unit vectors spaning the
% cone.
phi = zeros(nrClouds,1);

for k = 1:length(pointClouds)
    idx = getOcclusionPoints2D(pointClouds{k});
    
    % Get the cone-spaning vectors and turn them into unit vectors
    u1 = pointClouds{k}(idx(1),:);
    u1 = u1./norm(u1,2);
    
    u2 = pointClouds{k}(idx(2),:);
    u2 = u2./norm(u2,2);
    
    uv1(k,:) = u1;
    uv2(k,:) = u2;
    
    % Angle between the two unit vectors
    phi(k) = acos(dot(u1,u2));
    
    % Save unit vector to target
    m = mean(pointClouds{k});
    tgVector(k,:) = m./norm(m,2);
    
    % distance to target + margin
    distance(k) = norm(m, 2) + 0.5;    
end

% Allocate space for the detectionProbabilities
Pd = ones(nrClouds,1);

for k = 1:nrClouds
    % Make duplicat rows of the target unit vector
    tg = repmat(tgVector(k,:),nrClouds,1); % Extend matrix
    
    % Calculate the angle between the target vector and each cone-spanning
    % unit vector.
    angle1 = acos(dot(tg, uv1, 2));
    angle2 = acos(dot(tg, uv2, 2));
    
    totAngle = angle1 + angle2;
    
    % If a target is inside a cone, then the sum of the two angles between
    % the target and each unit vector will be the same as the total angle
    % between the cone-spanning unit vectors.
    tol = 1e-10;
    bounds = abs(totAngle - phi) < tol;
    
    idx = bounds & (distance(k) > distance);
    
    % Occlusion by several vehicles
    if sum(idx) > 1
        % Find occluding vehicle closest to EGO
        v = min(distance(idx));
        idx  = (distance == v);
    end
    
    if sum(idx) == 1
        % Take out the active cone-spanning unit vectors
        validVector1 = uv1(idx,:);
        validVector2 = uv2(idx,:);
        
        % distance between the cone-spanning unit vectors
        d = norm(validVector1 - validVector2,2);
        
        % The distance between the target unit vector to each of the cone
        % spanning unit vectors, use the smallest one to determine
        % detection probability
        d1 = norm(validVector1 - tgVector(k,:),2);
        d2 = norm(validVector2 - tgVector(k,:),2);
        
        dOrdered = sort([d1 d2]);
        
        % Normalize the smallest distance to a scale of 0-10
        % A symmetric exponential function is used to model the probability
        % of detection, with a lower bound of L. Alpha determines how fast
        % the detection probability decreases when moving to the middle of
        % the cone.
        x = mat2gray(dOrdered(1),[0 d])*10;
        a = 0;
        b = 10;
        alpha = 2;
        L = 0.1;
        f = @(x) L + (1-L)*(exp(-alpha*(x-a)) + exp(alpha*(x-b)));
        Pd(k) = f(x);
        
        if Pd(k) > 1
            Pd(k) = 1;
        end
    end
end
end