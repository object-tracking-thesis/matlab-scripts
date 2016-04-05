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

%% print it to check that everything went correctly
frame = [staticFrame(:,1:3); lidarData{1}(:,1:3)];
minx = min(frame(:,1))-5;
miny = min(frame(:,2))-5;
minz = min(frame(:,3))-5;

maxx = max(frame(:,1))+5;
maxy = max(frame(:,2))+5;
maxz = max(frame(:,3))+5;

staticCloud = pointCloud(frame(:,1:3));
player = pcplayer([minx maxx],[miny maxy], [minz maxz]);
view(player,staticCloud)  