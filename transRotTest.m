%% load lidar and oxts data
path1 = '../matlab-scripts/';
path1 = '~/Google Drive/Thesis Work AF - Object Tracking in Point Cloud/Data/2011_09_26_drive_crossing_dynamic/';
path2 = 'oxts/data/';
path3 = 'kitti/';
path4 = 'oxts/dataformat.txt';

oxtsData = loadOxtsDir(strcat(path1,path2));
lidarData = loadLidarDir(strcat(path1,path3));
oxtsStruct = oxts2struct(oxtsData,strcat(path1,path4));
Num = length(lidarData);

%% create the static frame
s = 1; %what frame to use as static?
staticFrame = lidarData{s}(:,1:3);
staticGPS = [oxtsStruct.lat(s) oxtsStruct.lon(s)];
staticAlt = oxtsStruct.alt(s);
staticYaw = oxtsStruct.yaw(s);

%% translate the live frames to the static environment
liveFrames = cell(1,Num);
for i=1:Num
    liveGPS = oxtsData{i}(1:2);
    liveAlt = oxtsData{i}(3);
    liveFrames{i} = translateFrame(lidarData{i}(:,1:3), staticGPS, liveGPS, staticAlt, liveAlt);
end

%% rotate the live frames to the static environment
for i=1:Num
    liveYaw = oxtsStruct.yaw(i);
    liveFrames{i} = rotateFrame(liveFrames{i},liveYaw,-staticYaw);
end

%% subtract the static part from all frames
cleanedFrames = cell(1,Num);
limit = 0.1;
for i=180:Num
    cleanedFrames{i} = subPC(staticFrame,liveFrames{i}(:,1:3),limit); 
end

%% Plot Animation 
wd = 50; % Axis size
   figure('Name','Removing static map',...
          'Position', [50 100 1200 500]);
for i = 180:Num
   
   subplot(1,2,1)
        lidarPlot(cleanedFrames{i})
        axis(wd*[-1 1 -1 1 -1/(0.5*wd) 1])
        grid off; box off;
        str = sprintf('Frame: %d / %d', [i length(lidarData)]);
        title(str);
   subplot(1,2,2)
        h = lidarPlot(lidarData{i});
        set(h,'MarkerEdgeColor', [1 0 0],'MarkerFaceColor', [1 0 0]);
        axis(wd*[-1 1 -1 1 -1/(0.5*wd) 1])
        grid off; box off;
        pause(0.25)
   
end

%% test plot difference between frame 200 and frame 180 translated to 200
frame200 = lidarData{201}(:,1:3);
frame180 = lidarData{181}(:,1:3);
frame180to200 = translateFrame(lidarData{181}(:,1:3), oxtsData{181}(1:2),...
                                oxtsData{201}(1:2), oxtsData{181}(3), oxtsData{201}(3));
figure;
lidarPlot(frame180,0.2)
hold on
%lidarPlot(frame180to200,0.5)
hold on
lidarPlot(frame200,0.8)
axis(wd*[-1 1 -1 1 -1/(0.5*wd) 1])
grid off; box off;

%% test diff on the forward and leftwards velocity
velx = [];
vely = [];
for i=180:200
    velx = [velx oxtsData{i}(12)];
    vely = [vely oxtsData{i}(13)];
end
[sum(velx)/10 sum(vely)/10]