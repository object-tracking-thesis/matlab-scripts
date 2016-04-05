%% clear
clear all;
clc;

%% load lidar and gps data for live
path1 = '~/Downloads/thesis/share/pcap_scenarios/';
%path2 = '20150204_WKS128_CYC2/';
%path2 = '20150204_WKS128_CYC3/';
%path2 = '20150204_WKS128_CYC9/';
%path2 = '20150204_WKS128_CYC11/';
%path2 = '20150204_WKS128_PED6/';
%path2 = '20150204_WKS128_PED8/';
path2 = 'car/';
path3 = 'gps/';
path4 = 'pcd/';

[affineStart, velStart, accStart, affineEnd, velEnd, accEnd, offset] = ... 
    loadGPSDir(strcat(path1,path2,path3));
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

%% downsample the static map
staticCloud = pointCloud(staticFrame(:,1:3));
gridStep = 0.8;
staticCloudDownsampled = pcdownsample(staticCloud, 'gridAverage', gridStep);

%% transform the live frames via their affine matrix and remove the offset to static
liveFrames = cell(1,Num);
for i=1:Num
    liveFrames{i} = transformFrameTransMat(lidarData{i}(:,1:3), affineStart{i});
    angle = deg2rad(4);
    rotmat = [cos(angle) sin(angle) 0; -sin(angle) cos(angle) 0; 0 0 1];
    newOffset = -offset{i};%+[-9 4 0]';
    liveFrames{i} = transformFrameTransMat(liveFrames{i}(:,1:3), [rotmat newOffset;zeros(1,3) 1]);
end

%% show live frames
figure
for frame=1:Num
    plot3(0,0,0,'.','MarkerSize',20);
    hold on
    live = pointCloud(liveFrames{frame}(:,1:3));
    pcshow(live)
    hold off
    str = sprintf('Frame: %d', [frame]);
    title(str);
    zoom(2)
    pause(0.2)
end

%% plot difference between transformed and untransformed liveframes
figure
for frame=1:Num
    plot3(0,0,0,'.','MarkerSize',20);
    hold on
    live = pointCloud(lidarData{frame}(:,1:3));
    liveTrans = pointCloud(liveFrames{frame}(:,1:3));
    pcshowpair(live, liveTrans)
    hold off
    zoom(2)
    pause(0.5)
end

%% plot difference between static and live map
figure
for frame=1:Num
    plot3(0,0,0,'.','MarkerSize',20);
    hold on
    liveCloud = pointCloud(liveFrames{frame}(:,1:3));
    pcshowpair(staticCloudDownsampled, liveCloud)
    hold off
    zoom(4)
    pause(0.3)
end

%% test plots
figure
live = pointCloud(lidarData{1}(:,1:3));
angle = pi;
rotmat = [cos(angle) sin(angle) 0; -sin(angle) cos(angle) 0; 0 0 1];
liveTrans = pointCloud(transformFrameTransMat(lidarData{1}(:,1:3),[rotmat [100 0 0]'; zeros(1,3) 1]));
pcshowpair(live, liveTrans)
zoom(2)