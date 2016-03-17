%% test MHTFinstance with twoTracks.m

run('twoTracks.m');
%%
% Setup measurements for two of the targets, x1_m & x2_m, add in no
% clutter. NO OPTIMAL ASSIGNMENT ALGO USED 

T1m = x1_m(1:2,:);
T2m = x2_m(1:2,:);

rawscan = cell(1,N);

for j = 1:N
    rawscan{j} = [T1m(:,j) T2m(:,j) clutterMeas{j}]; % clutterMeas{j}(:,1:randi(5))
    rawscan{j} = rawscan{j}(:, randperm(length(rawscan{j}))); % Shuffle measurement order (columns)
end

%% Place measurements in scanobject vector
clear Scans
Scans(1,N) = Scan;

for j = 1:N
    Scans(j).addMeasurements(rawscan{j});
end


%% Create MHTF instance
clear bestHypos
clear bestTracks
nrHypos = 15;
bestHypos(1,N) = Hypothesis;
bestTracks(1,N) = Tracks;

tic
for k = 1:N
    disp(k)    
    if k == 1        
        MHTF = MHTFinstance(nrHypos,N,Scans(k));
        bestHypos(k) = MHTF.bestHypo.copy();
        bestTracks(k) = bestHypos(k).tracks.copy();
        
        if isempty(bestTracks(k).trackId)
            bestTracks(k) = Tracks;
        else
            for j = 1:length(bestTracks(k).trackId)
                bestTracks(k).track(j) = bestTracks(k).track(j).copy();
            end
        end

    else
        % Now we iterate through each hypothesis and investigate
        MHTF.iterate(Scans(k));
        bestHypos(k) = MHTF.bestHypo.copy();
        bestTracks(k) = bestHypos(k).tracks.copy();
        
        if isempty(bestTracks(k).trackId)
            bestTracks(k) = Tracks;
        else
            for j = 1:length(bestTracks(k).trackId)
                bestTracks(k).track(j) = bestTracks(k).track(j).copy();
            end
        end
        
    end
end
toc


%% Look at best hypo history 
if 1
Targets = cell(1,N);

for k = 1:N
    try
    Targets{k} = [bestTracks(k).track.expectedValue];
    catch e
       Targets{k} = [bestTracks(k).track]; 
    end

end

figure(h1)
filename = 'testnew51.gif';
for k = 1:N
    hold on
    %p3 = plot(T1(1,k),T1(3,k),'sk','MarkerFaceColor','k');
    try
    plot(Targets{k}(1,:),Targets{k}(3,:),'sk','MarkerFaceColor','k')
    catch e
    end
    %lg = legend([p1 p11 p2 p21 p3],'Target 1','Target 1 measurement', 'Target 2', 'Target 2 measurement','MHTF state estimates');
    %lg.FontSize = 14;
    hold off    
    txt = sprintf('t = %d (1:%d)',[k N]);
    title(txt,'FontSize',14,'FontWeight','Bold')
    
%               drawnow
%       frame = getframe(1);
%       im = frame2im(frame);
%       [imind,cm] = rgb2ind(im,256);
%       if k == 1;
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf)
%       else
%           imwrite(imind,cm,filename,'gif','WriteMode','append')
%       end
    pause(0.1)
    
end
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

