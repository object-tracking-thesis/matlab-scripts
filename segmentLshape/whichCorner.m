% prototype for a function to return which corner of the car we are
% actually looking at 
%
%       ^ 
% UL ___|___ UR
%   |   |   |
%   |       |
%   |       |
%   |       |
%   |       |
% LL|______ |LR
%

% We are at time k, only get x & y
function corner = whichCorner(egoPosNow, egoPosPrev, egoRotNow, egoRotPrev, clusterNow, clusterPrev)


R = @(phi)[cos(phi) -sin(phi); sin(phi) cos(phi)];

% Compensate for rotations on egoPosition 
egoPosPrev = egoPosPrev' * R(egoRotPrev); % Previous position now in earth frame
egoPosNow = egoPosNow' * R(egoRotNow); % Current position now in earth frame 

% Calculate position diff 

posDiff = egoPosNow - egoPosPrev; % Difference in ego in earthfame  

disp(posDiff)

% Compensate target movements 

mean(clusterNow*R(egoRotNow),2)'
mean(clusterPrev*R(egoRotPrev),2)'


targetDiff = mean(clusterNow*R(egoRotNow),2)' + posDiff - mean(clusterPrev*R(egoRotPrev),2)';

corner = 1;



negComponent = ['UL', 'UR'];
posComponent = ['LL', 'LR'];


 



end