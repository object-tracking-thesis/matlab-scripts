%% test MHTFinstance with twoTracks.m

run('twoTracks.m');
%%
% Setup measurements for two of the targets, x1_m & x2_m, add in no
% clutter. NO OPTIMAL ASSIGNMENT ALGO USED 

T1m = x1_m(1:2,:);
T2m = x2_m(1:2,:);

rawscan = cell(1,N);

for j = 1:N
    rawscan{j} = [T1m(:,j) T2m(:,j)]; % clutterMeas{j}(:,1:randi(5))
    rawscan{j} = rawscan{j}(:, randperm(length(rawscan{j}))); % Shuffle measurement order (columns)
end

%% Place measurements in scanobject vector
clear Scans
Scans(1,N) = Scan;

for j = 1:N
    Scans(j).addMeasurements(rawscan{j});
end


%% Create MHTF instance
nrHypos = 3;
bestHypos(1,N) = Hypothesis;
bestTracks(1,N) = Tracks;
for k = 1:N
    if k == 1
        % We have 2 hypos that have 2 targets, and one hypo with 1 target
        MHTF = MHTFinstance(nrHypos,5,Scans(k));
        bestHypos(k) = MHTF.bestHypo.copy();
        bestTracks(k) = bestHypos(k).tracks.copy();
        
        for j = 1:bestTracks(k).trackId(end)
           bestTracks(k).track(j) = bestTracks(k).track(j).copy();
        end
    else
        % Now we iterate through each hypothesis and investigate
        MHTF.iterate(Scans(k));
        bestHypos(k) = MHTF.bestHypo.copy();
        bestTracks(k) = bestHypos(k).tracks.copy();
        
        for j = 1:bestTracks(k).trackId(end)
            bestTracks(k).track(j) = bestTracks(k).track(j).copy();
        end
        
    end
    fprintf('k = %d, length(hypoStorage) = %d, length(tempStorage) = %d\n', [k, length(MHTF.hypoStorage),length(MHTF.tempStorage)]);
    fprintf('')
end



%% Look at best hypo history 
T1 = ones(4,N);
T2 = ones(4,N);


for k = 1:N
    T1(:,k) = bestTracks(k).track(1).expectedValue;
    T2(:,k) = bestTracks(k).track(2).expectedValue;

end

for k = 1:N
    figure(h1)
    hold on
    plot(T1(1,k),T1(3,k),'sk','MarkerFaceColor','k')
    plot(T2(1,k),T2(3,k),'sk','MarkerFaceColor','k')
    hold off
    pause(2)
end




%% DEBUG FUNCTIONS
DEBUG = 0;
if DEBUG
    k = 10;
    fprintf('MEASUREMENTS AT K = %d\n',k)
    disp([T1m(:,k), T2m(:,k) ])
    fprintf('STATE AT K = %d\n',k)
    disp([x1_state(:,k), x2_state(:,k) ])
    fprintf('-------------------------------\n')
    % DEBUG PRINTER FILTER
    
    fprintf('MHTF.hypoStorage(1) \n\n')
    disp([MHTF.hypoStorage(1).tracks.track.expectedValue])
    fprintf('MHTF.hypoStorage(2) \n\n')
    disp([MHTF.hypoStorage(2).tracks.track.expectedValue])
    fprintf('MHTF.hypoStorage(3) \n\n')
    disp([MHTF.hypoStorage(3).tracks.track.expectedValue])
    % TEMPSTORAGE 1:3
    fprintf('MHTF.tempStorage(1) \n\n')
    disp([MHTF.tempStorage(1).tracks.track.expectedValue])
    fprintf('MHTF.tempStorage(2) \n\n')
    disp([MHTF.tempStorage(2).tracks.track.expectedValue])
    fprintf('MHTF.tempStorage(3) \n\n')
    disp([MHTF.tempStorage(3).tracks.track.expectedValue])
    % TEMPSTORAGE 4:6
    fprintf('MHTF.tempStorage(4) \n\n')
    disp([MHTF.tempStorage(4).tracks.track.expectedValue])
    fprintf('MHTF.tempStorage(5) \n\n')
    disp([MHTF.tempStorage(5).tracks.track.expectedValue])
    fprintf('MHTF.tempStorage(6) \n\n')
    disp([MHTF.tempStorage(6).tracks.track.expectedValue])
    % TEMPSTORAGE 7:9
    fprintf('MHTF.tempStorage(7) \n\n')
    disp([MHTF.tempStorage(7).tracks.track.expectedValue])
    fprintf('MHTF.tempStorage(8) \n\n')
    disp([MHTF.tempStorage(8).tracks.track.expectedValue])
    fprintf('MHTF.tempStorage(9) \n\n')
    disp([MHTF.tempStorage(9).tracks.track.expectedValue])
end

