%% MHT filter implementation for twoTracks.m
% MAKE SURE TO RUN twoTracks.m FIRST! 
% No clustering, pruning or gating. No new tracks, no dead tracks, no
% clutter.

run('twoTracks.m')

%% Movement and Measurement model knowledge

T = 1;
A = [1 T 0 0;
     0 1 0 0;
     0 0 1 T;
     0 0 0 1];

Q_mu = [0 0 0 0]';

Q_sigma = 0.05.*[0 0 0 0;
                 0 1 0 0;
                 0 0 0 0;
                 0 0 0 1];

H = [1 0 0 0;
     0 0 1 0];

R_mu = [0 0]';

R_sigma = 0.125.*[1 0;
                  0 1];

%% MHT loop
% These values are set when twoTracks2.m is run
        % Set init values for our two tracks (we know we have two tracks & where they are)
        % N = 5;

        % x1_init = [1 1]';
        % x2_init = [1 5]';

%% Set up Measurement Sequence

% Containes the measurements gathered for each scan (N scans)
Z = ones(2,2,N);

for k = 1:N
    Z(:,:,k) = [x1_m(:,k) x2_m(:,k)];
end
%% Set up init state vector for predictions

s_mu = zeros(4,2);
s_P = zeros(4,4,2);
%% Struct storage for filter componentes in each time step
storage = cell(1,N); % storage for generated hypotheses

%% Run loop

for k = 1:N
    if k == 1
        % Predict moments (mu & P)
        s_mu(:,1) = A*x1_init;
        s_P(:,:,1) = A*zeros(4,4)*A' + Q_sigma;
        
        s_mu(:,2) = A*x2_init;
        s_P(:,:,2) = A*zeros(4,4)*A' + Q_sigma;
        
        % Formulate hypotheses, should be changed to combvec
        hypMat = perms(fliplr(1:length(Z(:,:,k))));
        
        [r,c] = size(Z(:,:,k)); % Number of scans and measurements in each
        hypNum = factorial(r); % Number of hypotheses, should be changed to use combvec
        
        G = cell(hypNum,1);
        % Calculate Posterior Moments for each Hypothesis
        for hyp = 1:hypNum % Choses active hypothesis
            h = hypoStorage; 
            h.predMu = s_mu;
            h.predCov = s_P;
            h.hypo = hypMat(hyp,:);
            h.beta = 1;
            h.postCov = ones(4,4,c);
            h.postMu = ones(4,1,c);
            h.parentHypo = 'root';
            for tg = 1:c % Choses target to calc post. for
               
                h.postMu(:,1,tg) = h.predMu(:,tg) + h.predCov(:,:,tg)*H'*inv(H*h.predCov(:,:,tg)*H' + R_sigma)*(Z(:,h.hypo(tg),k) - H*h.predMu(:,tg));
                h.postCov(:,:,tg) = h.predCov(:,:,tg) - h.predCov(:,:,tg)*H'*inv(H*h.predCov(:,:,tg)*H' + R_sigma)*H*h.predCov(:,:,tg);
               
                % Calculate Beta
                h.beta = h.beta*mvnpdf(Z(:,h.hypo(tg),k), H*h.predMu(:,tg), H*h.predCov(:,:,tg)*H' + R_sigma);
                
            end
            G{hyp,1} = h;
        end
        
        % Calculate sum of beta for k=1 hypotheses 
        totalBeta = 0;
        for j=1:length(G)
            totalBeta = totalBeta + G{j}.beta;
        end
        
        % Calculate and set alpha
        for j=1:length(G)
            G{j}.alpha = G{j}.beta/totalBeta;
        end
        
        % Store components 
        storage{k} = G;
        
    else
        
        rootSumBeta = 0; % Used to sum total beta for the branches of the root
        % select root hypotheses
        rootHypos = storage{k-1};
        % Nr of hypotheses we have at root, start out with first
        for j = 1:length(rootHypos)
            s_mu = ones(4,1,2);
            s_P = ones(4,4,2);
            % Predict moments (mu & P)
            s_mu(:,1,1) = A*rootHypos{j}.postMu(:,1,1);
            s_P(:,:,1) = A*rootHypos{j}.postCov(:,:,1)*A' + Q_sigma;
            
            s_mu(:,1,2) = A*rootHypos{j}.postMu(:,1,2);
            s_P(:,:,2) = A*rootHypos{j}.postCov(:,:,2)*A' + Q_sigma;
            
            % Formulate new hypotheses
            hypMat = perms(fliplr(1:length(Z(:,:,k))));
            [r,c] = size(Z(:,:,k)); % Number of scans and measurements in each
            hypNum = factorial(r); % Number of hypotheses, should be changed to combvec
            
            G = cell(hypNum,1);
            % Calculate Posterior Moments for each Hypothesis
            for hyp = 1:hypNum % Choses active hypothesis
                h = hypoStorage;
                h.predMu = s_mu;
                h.predCov = s_P;
                h.hypo = hypMat(hyp,:);
                h.beta = 1;
                h.postCov = ones(4,4,c);
                h.postMu = ones(4,1,c);
                h.parentHypo = rootHypos{j}.hypoNr;
                for tg = 1:c % Choses target to calc post. for

                    h.postMu(:,1,tg) = h.predMu(:,tg) + h.predCov(:,:,tg)*H'*inv(H*h.predCov(:,:,tg)*H' + R_sigma)*(Z(:,h.hypo(tg),k) - H*h.predMu(:,tg));
                    h.postCov(:,:,tg) = h.predCov(:,:,tg) - h.predCov(:,:,tg)*H'*inv(H*h.predCov(:,:,tg)*H' + R_sigma)*H*h.predCov(:,:,tg);
                    
                    % Calculate Beta
                    h.beta = h.beta*mvnpdf(Z(:,h.hypo(tg),k), H*h.predMu(:,tg), H*h.predCov(:,:,tg)*H' + R_sigma);
                end
                h.beta = h.beta * rootHypos{j}.alpha;
                rootSumBeta = rootSumBeta + h.beta;
                G{hyp,1} = h;
            end
            % Place G in storage 
            storage{k} = [storage{k}; G];
            
        end
        % Calculate and set alpha
        for n=1:length(storage{k})
            storage{k}{n}.alpha = storage{k}{n}.beta/rootSumBeta;
        end
    end
end


%% store alpha probability for each hypo 
alphas = cell(1,N);
for j = 1:length(storage)
    for k = 1:length(storage{j})
        alphas{1,j} = [alphas{1,j};storage{j}{k}.alpha];
    end
end

% find max probable alpha 
idx = ones(1,length(alphas));
for j = 1:length(alphas)
    [r, c] = max(alphas{j});
    idx(j) = c;
end

% Get best hypo posterior mu

bestPostT1 = ones(4,1,N);
bestPostT2 = ones(4,1,N);

for k = 1:length(storage)
   bestPostT1(:,:,k) = storage{k}{idx(k)}.postMu(:,:,1);
   bestPostT2(:,:,k) = storage{k}{idx(k)}.postMu(:,:,2);
    
end

    
%% Plot most likely hypothesis

figure(h1)
hold on

plot(bestPostT1(1,:),bestPostT1(3,:),'db','MarkerFaceColor','b')
plot(bestPostT2(1,:),bestPostT2(3,:),'dr','MarkerFaceColor','r')

