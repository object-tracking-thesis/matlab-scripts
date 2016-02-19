path = 'static/wsp_square_3.txt';
path2 = 'static/wsp_square_4.txt';
staticFrame1 = xyzirgb2matlab(path);
staticFrame2 = xyzirgb2matlab(path2);
staticFrame = [staticFrame1; staticFrame2];


%%
minx = min(staticFrame(:,1))-5;
miny = min(staticFrame(:,2))-5;
minz = min(staticFrame(:,3))-5;

maxx = max(staticFrame(:,1))+5;
maxy = max(staticFrame(:,2))+5;
maxz = max(staticFrame(:,3))+5;

staticCloud = pointCloud(staticFrame(:,1:3));
player = pcplayer([minx maxx],[miny maxy], [minz maxz]);
view(player,staticCloud)  