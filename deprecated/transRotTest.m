%% load lidar and oxts data
path1 = '../matlab-scripts/';
path2 = 'oxts/';
path3 = 'pcd/';
path4 = 'dataformat.txt';

oxtsData = loadOxtsDir(strcat(path1,path2));
lidarData = loadLidarDir(strcat(path1,path3));
oxtsStruct = oxts2struct(oxtsData,strcat(path1,path2,path4));
Num = length(lidarData);

%% create the static frame
s = 1; %what frame to use as static?
staticFrame = lidarData{s}(:,1:3);
staticGPS = [oxtsStruct.lat(s) oxtsStruct.lon(s)];
staticAlt = oxtsStruct.alt(s);
staticYaw = oxtsStruct.yaw(s);

%% rotate the live frames to the static environment
liveFrames = cell(1,Num);
for i=1:Num
    liveYaw = oxtsStruct.yaw(i);
    liveFrames{i} = rotateFrame(lidarData{i}(:,1:3),liveYaw,-staticYaw);
end

%% translate the live frames to the static environment
for i=1:Num
    liveGPS = oxtsData{i}(1:2);
    liveAlt = oxtsData{i}(3);
    liveFrames{i} = translateFrame(liveFrames{i}, staticGPS, liveGPS, staticAlt, liveAlt);
end


%% subtract the static part from all frames
cleanedFrames = cell(1,Num);
limit = 0.5;
for i=180:Num
    cleanedFrames{i} = subPC(staticFrame,liveFrames{i}(:,1:3),limit); 
end

%% Plot Animation 
wd = 50; % Axis size
   h = figure('Name','Removing static map',...
          'Position', [50 100 1200 500]);
   h.Renderer = 'opengl';
for i = 180:Num
   
   subplot(1,2,1)
        lidarPlot(cleanedFrames{i})
        %hold on
        %scatter3(0,0,0,'MarkerEdgeColor',[0 1 0],'MarkerFaceColor',[0 1 0])
        %hold off
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

%% test plot difference between frame 180 and frame X translated to 180
x = 200;
frameX = lidarData{x}(:,1:3);
frame180 = lidarData{181}(:,1:3);
frameXto180 = translateFrame(lidarData{x}(:,1:3), oxtsData{181}(1:2),...
                                oxtsData{x}(1:2), oxtsData{181}(3), oxtsData{x}(3));
frameXto180 = rotateFrame(frameXto180,oxtsStruct.yaw(x),-oxtsStruct.yaw(181));
%limit = 1;
%frameXto180 = subPC(frame180,frameXto180,limit);

figure;
lidarPlot(frame180,[0.5 0 0])
hold on
lidarPlot(frameXto180,[0 0.5 0])
hold on
lidarPlot(frameX,[0 0 0.5])
legend('180','Xto180','X')
axis(wd*[-1 1 -1 1 -1/(0.5*wd) 1])
grid off; box off;

%% test plot conjunction between before crossing and after
wd = 100;
frame180 = cleanedFrames{181}(:,1:3);
frame269 = cleanedFrames{210}(:,1:3);

figure;
lidarPlot(frame180,0.2)
hold on
lidarPlot(frame269,0.8)
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