   
%% Monte Carlo Simulation of one transition through f. 
% Do this to get initial covariance vector.
N = 10000;

X = mvnrnd(zeros(7,1), eye(7), N);
Y = zeros(size(X));

for j = 1:N
    Y(j,:) = f(X(j,:))'; % From gaussian to first state
    for h = 1:1; % The non-gaussian to second nongaussian state
        Y(j,:) = f(Y(j,:))';
    end
end

% 

P = cov(Y);
eig(P)