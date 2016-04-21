%% Solve optimal L-shape
% S - set of all points 
% u = (c1,c2,n1,2), Lp = n1*x + n2*y + c1 = 0 & Lq = -n2*x + n1*y + c2 = 0

rng(1337)
% Generate som test data 
m = 100; % nr of total points 
r = 1:m/2;
alongX = mvnrnd([1 0], [0.2 0.01;0.01 0.01],m/2); %[r', zeros(m/2,1)];%
alongY = mvnrnd([0 1.5], [0.01 0.01;0.01 0.4],m/2); %[zeros(m/2,1), r'];%

% Sort data to be ordered counterclockwise 
alongX = flip(sortrows(alongX,1));
alongY = sortrows(alongY,2);

S = [alongX; alongY];



phi = pi/6;

R = [cos(phi) -sin(phi); sin(phi) cos(phi)];   

S = S*R;

%S = S - repmat([10 1],m,1);

f = figure;
    hold on; grid on; box on;
    p1 = plot(S(1:0.5*m,1), S(1:0.5*m,2), 'rx');
    p2 = plot(S(0.5*m+1:end,1), S(0.5*m+1:end,2), 'bx');
    plot(0,0,'sk')
    ax = 10;
    axis([-ax ax -ax ax])
    legend([p1 p2], 'Along Y', 'Along X')
%%
if 0
f = figure;

for k = 1:length(S)    
    clf
    figure(f)
    hold on
    plot(S(:,1), S(:,2), 'bx')
    plot(S(k,1), S(k,2),'rx')
    
    pause(1)
    axis equal
end
    
end
%% This is where the algorithm starts
tic
v = 1; % initial index 
m = size(S,1); % Nr of rows/data
lambdaOp = 1e10; % approx inf, Op stands for Optimal 

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
%%
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
    v
    [n, lambda] = eig(M22 - M12'*inv(M11)*M12); % M22 - M12'*inv(M11)*M12
    % Keep smallest eigenvalue & corresponding vetor 
    [d, i] = min(diag(lambda));
    lambda = d;
    n = n(:,i);
    
    if lambdaOp > lambda   
       disp(v)
       c = -1.*inv(M11)*M12*n;
       lambdaOp = lambda; 
       uOp = [c;n];
    end
    v = v + 1;    
end
toc
% Done
%% Draw lines from uOp

c1 = uOp(1);
c2 = uOp(2);
n1 = uOp(3);
n2 = uOp(4);

% Define the two perp lines
x = -ax:ax;
y1 = -n1/n2*x - c1/n2;
y2 = n2/n1*x - c2/n1;

figure(f)

plot(x,y1,'-g')
plot(x,y2,'-g')

xc = -n1*c1 + n2*c2;
yc = -n2*c1 - n1*c2;
plot(xc, yc, 'sk','MarkerFace','k')

% 

vec1 = [xc yc;1 n2/n1-c2/n1];
vec2 = [xc yc;1 -n1/n2-c1/n2];

plot(vec1(:,1), vec1(:,2),'--k')
plot(vec2(:,1), vec2(:,2),'--m')


%% Prob distro 

% [H, prH] = ISEDVPE(S, 1);

% plot(H(2,:),prH,'x')










