% ISED ("Incremental Sub-matrix Eigen Decomposition") method
% Takes a set of xy points and returns coefficients of the optimal
% orthogonal lines that span an optimal L-share over the data (in least 
% square terms). NOTE! In order for method to work xy points need to be 
% ordered clockwise/counterclockwise. The lines are defined as: 
%
% Yp(x) = -(n1*x + c1)/n2  
% Yq(x) =  (n2*x - c2)/n1
% 
% Function:
%        uOp = ISED(S)
%   [uOp, i] = ISED(S)
%
% Parameters:
%     S - mx2 vectors with xy points. 
%   uOp - vector containing [c1 c2 n1 n2]
%     i - index in vector S where the orthogonal split occurs 

function [uOp, varargout] = ISED(S)

v = 1; % initial index 
m = size(S,1); % Nr of rows/data
lambdaOp = 1e10; % approx inf, Op stands for Optimal 
idx = 1; % Orthogonal split index 

% calc M
p = v;
q = m-v;

Xp = S(1:v,1); % Xp
Yp = S(1:v,2); % Yp

Yq = S(v+1:end,2); % Yq
Xq = S(v+1:end,1); % Xq

b1 = sum(Xp); % sumXp
c1 = sum(Yp); % sumYp

b2 = sum(Yq); % sumYq 
c2 = -1*sum(Xq); % sumXq

d3 = sum(Xp.^2) + sum(Yq.^2); % sqSumXpYq
c3 = sum(Xp.*Yp) - sum(Xq.*Yq); % sqDiffXpYq
d4 = sum(Yp.^2) + sum(Xq.^2); % sqSumYpXq 

M = [ p  0 b1 c1;...
      0  q b2 c2;...
     b1 b2 d3 c3;... 
     c1 c2 c3 d4];
 
M12 = M(1:2,3:4);
M11 = M(1:2,1:2);
M22 = M(3:4,3:4);

[n, lambda] = eig(M22 - M12'*(M11\M12)); % M22 - M12'*inv(M11)*M12
% Keep smallest eigenvalue & corresponding vetor 
[d, i] = min(diag(lambda));
lambda = d;
n = n(:,i);

if lambdaOp > lambda
   c = -1.*inv(M11)*M12*n;
   lambdaOp = lambda; 
   uOp = [c;n];
end

while v < (m-1)
    % calculate deltaM 
    b1 = S(v+1,1); % x_{v+1}
    c1 = S(v+1,2); % y_{v+1}
    b2 = -c1; 
    c2 = b1;
    d3 = b1^2 - c1^2; % x_{v+1}^{2} - y_{v+1}^{2}
    c3 = 2*b1*c1;     % 2 * x_{v+1} * y_{v+1}
    d4 = c1^2 - b1^2; % y_{v+1}^{2} - x_{v+1}^{2}

    dM = [1  0  b1 c1;...
          0 -1  b2 c2;...
          b1 b2 d3 c3;...
          c1 c2 c3 d4];
      
    M = M + dM;
    
    M12 = M(1:2,3:4);
    M11 = M(1:2,1:2);
    M22 = M(3:4,3:4);
    
    [n, lambda] = eig(M22 - M12'*(M11\M12)); % M22 - M12'*inv(M11)*M12
    % Keep smallest eigenvalue & corresponding vetor 
    [d, i] = min(diag(lambda));    
    lambda = d;
    n = n(:,i);
    
    if lambdaOp > lambda      
       idx = v;       
       c = -1.*inv(M11)*M12*n;
       lambdaOp = lambda; 
       uOp = [c;n];
    end
    v = v + 1;    
end

if nargout == 2     
    varargout{1} = idx;
end


end