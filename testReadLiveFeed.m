%% clear
clear all;
clc;

%% load lidar and gps data for live
path1 = '~/Downloads/thesis/share/pcap_scenarios/';
path2 = '20150204_WKS128_CYC2/';
%path2 = '20150204_WKS128_CYC9/';
%path2 = '20150204_WKS128_PED6/';
%path2 = '20150204_WKS128_PED8/';
path3 = 'gps/';
path4 = 'pcd/';

affine = loadGPSDir(strcat(path1,path2,path3));
lidarData = loadLidarDir(strcat(path1,path2,path4));
Num = length(lidarData);

%% load lidar and gps data for static
path1 = 'static/wsp_square_1.txt';
path2 = 'static/wsp_square_2.txt';
path3 = 'static/wsp_square_3.txt';
path4 = 'static/wsp_square_4.txt';
staticFrame1 = xyzirgb2matlab(path1);
staticFrame2 = xyzirgb2matlab(path2);
staticFrame3 = xyzirgb2matlab(path3);
staticFrame4 = xyzirgb2matlab(path4);
staticFrame = [staticFrame1; staticFrame2; staticFrame3; staticFrame4];
clear staticFrame1;
clear staticFrame2;
clear staticFrame3;
clear staticFrame4;
staticGPS = load('static/staticgps.mat');
%staticOffset = [15.5 103.46 0]; %read visually from a plot
%staticGPS.staticGPS = staticGPS.staticGPS - staticOffset;

%% rotate and translate the live Frames to align with the static frame
liveFrames = cell(1,Num);
for i=1:Num
    liveFrames{i} = transformFrameTransMat(lidarData{i}(:,1:3), affine{i}, staticGPS.staticGPS);
end

%% subtract the static part from all frames
cleanedFrames = cell(1,Num);
limit = 0.5;
for i=250:250
    cleanedFrames{i} = subPC(staticFrame(:,1:3),liveFrames{i}(:,1:3),limit); 
end

%% print it to check that everything went correctly
frame = [staticFrame(:,1:3); liveFrames{140}(:,1:3)];
%frame = [cleanedFrames{150}(:,1:3)];
%frame = [cleanedFrames{250}(:,1:3)];
minx = min(frame(:,1))-5;
miny = min(frame(:,2))-5;
minz = min(frame(:,3))-5;

maxx = max(frame(:,1))+5;
maxy = max(frame(:,2))+5;
maxz = max(frame(:,3))+5;

staticCloud = pointCloud(frame(:,1:3));
player = pcplayer([minx maxx],[miny maxy], [minz maxz]);
view(player,staticCloud)