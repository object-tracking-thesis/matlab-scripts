%% load data
load thesisPlots/ped_ext_est.mat
load thesisPlots/ped_center_est.mat
load thesisPlots/group_ext_est.mat
load thesisPlots/group_center_est.mat


ext_est = group_ext_est;
center_est = group_center_est;
l = length(ext_est);
l = 29;
%% data prep
x_ext = [];
y_ext = [];
dx_ext = [];
dy_ext = [];
P_ext = zeros(3,3,l);
x_center = [];
y_center = [];
dx_center = [];
dy_center = [];
P_center = zeros(4,4,l);
for i= 1:l
    x_ext = [x_ext ext_est(i).state(1)];
    y_ext = [y_ext ext_est(i).state(2)];
    dx_ext = [dx_ext ext_est(i).state(3)];
    dy_ext = [dy_ext ext_est(i).state(4)];
    P_ext(:,:,i) = ext_est(i).P;
    x_center = [x_center center_est(i).mu(1)];
    y_center = [y_center center_est(i).mu(2)];
    dx_center = [dx_center center_est(i).mu(3)];
    dy_center = [dy_center center_est(i).mu(4)];
    P_center(:,:,i) = center_est(i).P;
end

%% plot positions
figure(1)
title('position estimates difference between extended targets and cluster center targets')
subplot(2,1,1);
plot(1:l, x_ext, '--b') 
hold on
plot(1:l, x_center, '--r')
xlabel('k')
ylabel('pos')
legend('x_{ext}', 'x_{center}')

state = 1;
sigEtt = 3*sqrt(reshape(P_ext(state,state,:), 1,l));
sigNo  = 3*sqrt(reshape(P_center(state,state,:), 1,l));

shadeX = [1:l, l:-1:1];
shadeY = [sigEtt+x_ext fliplr(-1.*sigEtt+x_ext)];
       
fl = fill(shadeX, shadeY, 'b');
fl.FaceAlpha = 0.05;
fl.EdgeColor = 'b';
fl.EdgeAlpha = 0;
       
shadeX = [1:l, l:-1:1];
shadeY = [sigNo+x_center fliplr(-1.*sigNo+x_center)];
       
fl = fill(shadeX, shadeY, 'r');
fl.FaceAlpha = 0.05;
fl.EdgeColor = 'r';
fl.EdgeAlpha = 0;

subplot(2,1,2)
plot(1:l, y_ext) 
hold on
plot(1:l, y_center) 
xlabel('k')
ylabel('pos')
legend('y_{ext}', 'y_{center}')

state = 2;
sigEtt = 3*sqrt(reshape(P_ext(state,state,:), 1,l));
sigNo  = 3*sqrt(reshape(P_center(state,state,:), 1,l));

shadeX = [1:l, l:-1:1];
shadeY = [sigEtt+y_ext fliplr(-1.*sigEtt+y_ext)];
       
fl = fill(shadeX, shadeY, 'b');
fl.FaceAlpha = 0.05;
fl.EdgeColor = 'b';
fl.EdgeAlpha = 0;
       
shadeX = [1:l, l:-1:1];
shadeY = [sigNo+y_center fliplr(-1.*sigNo+y_center)];
       
fl = fill(shadeX, shadeY, 'r');
fl.FaceAlpha = 0.05;
fl.EdgeColor = 'r';
fl.EdgeAlpha = 0;

%% plot velocities
figure(2)
title('velocity estimates difference between extended targets and cluster center targets')
subplot(2,1,1);
plot(1:l, dx_ext) 
hold on
plot(1:l, dx_center)
xlabel('k')
ylabel('vel')
legend('dx_{ext}', 'dx_{center}')

subplot(2,1,2)
plot(1:l, dy_ext) 
hold on
plot(1:l, dy_center) 
xlabel('k')
ylabel('vel')
legend('dy_{ext}', 'dy_{center}')