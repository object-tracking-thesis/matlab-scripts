% ISEDVPE ("Incremental Sub-matrix Eigen Decomposition Method for Vehicle 
% Pose Estimation") method.
% Takes a set of xy points and returns hypotheses with associated 
% probabilities over corner point (xc, yc) and orientation theta_i, where i
% = 0:3 since there are four possibilities of which direction the car is 
% going. 
%
% 
% Function:
%     [H, prH] = ISEDVPE(S, tau)
%
% Parameters:
%     S - mx2 vectors with xy points. 
%   tau - decay factor for prpbability of hypo (Pr = K.*exp(lambda/(tau*m^2)))
%     H - 6x(m-1) vector of hypotheses. Each colum is a hypo, where 
%         hypo = [xc yc theta_1 theta_2 theta_3 theta_4]';
%   prH - 1x(m-1) vector The associated probability of each hypothesis.

function [H, prH, varargout] = ISEDVPE(S,tau)

v = 1; % initial index 
m = size(S,1); % Nr of rows/data
uOp = zeros(4,m-1); % Line parameters 

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
n = n(:,i);
c = -1.*inv(M11)*M12*n;

H = zeros(6, m-1); % Each column is a hypothesis,
prH = zeros(1, m-1); % Each column is Prob of correspoinding hypo in H

xc = -n(1)*c(1) + n(2)*c(2);
yc = -n(2)*c(1) - n(1)*c(2);

theta = [1 1 1 1].*atan2(n(2),n(1)) + [0 1 2 3].*(pi/2);

h = [xc yc theta]';
pr = exp(-d/(tau*m^2));

K = pr; % Normalization factor 

H(:,1) = h;
prH(1) = pr;
uOp(:,1) = [c;n];

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
    v = v + 1;
    M12 = M(1:2,3:4);
    M11 = M(1:2,1:2);
    M22 = M(3:4,3:4);
    
    [n, lambda] = eig(M22 - M12'*(M11\M12)); % M22 - M12'*inv(M11)*M12
    % Keep smallest eigenvalue & corresponding vetor
    [d, i] = min(diag(lambda));
    n = n(:,i);
    c = -1.*inv(M11)*M12*n;
    
    xc = -n(1)*c(1) + n(2)*c(2);
    yc = -n(2)*c(1) - n(1)*c(2);
    
    theta = [1 1 1 1].*atan2(n(2),n(1)) + [0 1 2 3].*(pi/2);
    
    h = [xc yc theta]';
    pr = exp(-d/(tau*m^2));
    
    H(:,v) = h;
    prH(v) = pr;
    uOp(:,v) = [c;n];
    K = K + pr;

end

prH = prH./K; % Normalize to one

if nargout == 3    
    varargout{1} = uOp;
end

end