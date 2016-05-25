path = '/Users/markocotra/Google Drive/Thesis Work AF - Object Tracking in Point Cloud/Data/carClusters/carClusters_kitti_campus_01_50';
load(path);

allTg = carClusters;

tg1 = cell(1,1);
tg2 = cell(1,2);

for j = 1:length(allTg)

tg1{j} = allTg{j}{1};
disp(length(tg1{j}))
try
    tg2{j} = allTg{j}{2};
    disp(length(tg2{j}))
catch me
    tg2{j} = [];
    
end

end

%% 
fig = figure;
for j = 1:length(allTg)
%    plot(allTg{j}{1}(:,1), allTg{j}{1}(:,2), 'kx')
    if (j/2 - floor(j/2)) == 0 % even
        plot(tg1{j}(:,1), tg1{j}(:,2),'r.')
            hold on
            try
                plot(tg2{j}(:,1), tg2{j}(:,2), 'bx')
            catch me
            end
    else
        plot(tg1{j}(:,1), tg1{j}(:,2),'k.')
        hold on
        try
            plot(tg2{j}(:,1), tg2{j}(:,2), 'gx')
        catch me
        end
    end
    axis([0 40 -20 20])
  % axis equal
%     
     pause(0.3)
%     
    
    
end

%%
N = 40;

for j = 45
    hold off
    try
        plot3(tg1{j}(:,1), tg1{j}(:,2),tg1{j}(:,3), 'rx'); axis equal
    catch me
    end
    j
    pause(1)
end

%%
j = 23;
[m1, m2, uOp, filtNtg] = cornerPoint(tg2{j});
 
      c1 = uOp(1); c2 = uOp(2);
      n1 = uOp(3); n2 = uOp(4);
     
      xc = (-n1*c1 + n2*c2);
      yc = (-n2*c1 -n1*c2);
      plot(tg2{j}(:,1), tg2{j}(:,2), 'kx'); axis equal; hold on
      plot(filtNtg(:,1), filtNtg(:,2), 'r.'); axis equal; hold on
      
      plot(xc, yc,'cs','MarkerSize',20)
      
      for i = 1:4 % i = 2 is correct one for N = 1
          plot(4.*[0 cos(angle(i))]+[xc, xc], 4.*[0 sin(angle(i))]+[yc, yc], '-g','LineWidth',2)
          plot(0+xc, 0+yc, 'og','LineWidth',2)
      end
      
      
      
      angle = ones(1,4)*atan2(n2,n1) + [0 1 2 3]*pi/2;




