%% test MHTFinstance with twoTracks.m

run('twoTracks.m');
%%
% Setup measurements for two of the targets, x1_m & x2_m, add in no
% clutter. NO OPTIMAL ASSIGNMENT ALGO USED 

T1m = x1_m(1:2,:);
T2m = x2_m(1:2,:);

rawscan = cell(1,4);

for j = 1:4
    rawscan{j} = [T1m(:,j) T2m(:,j)]; % clutterMeas{j}(:,1:randi(5))
    rawscan{j} = rawscan{j}(:, randperm(length(rawscan{j}))); % Shuffle measurement order (columns)
end

% Place measurements in scanobject vector
Scans(1,4) = Scan;
for j = 1:4
    Scans(j).addMeasurements(rawscan{j});
end


%% Create MHTF instance
nrHypos = 3;
bestHypos(1,4) = Hypothesis;
for k = 1:4
    if k == 1
        % We have 2 hypos that have 2 targets, and one hypo with 1 target
        MHTF = MHTFinstance(nrHypos,1,Scans(k));
        bestHypos(k) = MHTF.bestHypo;
    else
        % Now we iterate through each hypothesis and investigate
        disp(k)
        MHTF.iterate(Scans(k));
        bestHypos(k) = MHTF.bestHypo;
        
    end
end

