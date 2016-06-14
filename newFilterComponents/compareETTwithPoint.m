%% Compare ETT with no ETT; 


A = load('/Users/markocotra//thesis/matlab-scripts/car1.mat');
data1 = A.car1adj;
data1 = data1(1:36);

% A = load('/Users/markocotra//thesis/matlab-scripts/car2.mat');
% data2 = A.car2adj;
%%
storage = struct('st', [],'Cov',[]);

car1ett = cell(1,length(data1));
car1no  = cell(1,length(data1));

car2ett = cell(1,length(car2));
car2no = cell(1,length(car2));

%% Car 1 ETT 

st1 = [14.75 0.75 13 6.1 0 1.8 4.5]'; % Car 1

ett = CarTarget;

ett.init(st1, [], [], 1);

fig = figure; fig.Position = [50 50 1600 800];
for k = 1:length(data1)

    ett.predict();
    
    ett.calcLikelihood(data1{k});
    
    [x, P] = ett.getState;
    
    car1ett{k} = struct('st', x,'Cov',P);
    
    hold off
    plot(data1{k}(:,1), data1{k}(:,2),'kx'); axis equal; hold on
    drawMyRide(x,'b');    
    drawMyRide(st1, 'c');
    Mxy = mean(data1{k}(:,1:2));
    
    title(sprintf('%i/%i',[k length(data1)]))
    
    plot(Mxy(1), Mxy(2),'gs','MarkerSize',10)
    waitforkey        

end

%% Car 2 (partial data) ETT 

st2 = [18.5 -13.4 4 pi/2 0 1.6 3.5]'; % Car 2

ett = CarTarget;

ett.init(st2, [], [], 1);

fig = figure; fig.Position = [50 50 1600 800];

% car2 is a mat file with complete sequence for car 2.
data2 = car2;

for k = 1:length(data2)

    ett.predict();
    
    ett.calcLikelihood(data2{k});
    
    [x, P] = ett.getState;
    
    car2ett{k} = struct('st', x,'Cov',P);
    
    hold off
    plot(data2{k}(:,1), data2{k}(:,2),'kx'); axis equal; hold on
    drawMyRide(x,'b');    
    %drawMyRide(st2, 'c');
    Mxy = mean(data2{k}(:,1:2));
    
    plot(Mxy(1), Mxy(2),'gs','MarkerSize',10)
    
    title(sprintf('%i/%i',[k length(data2)]))
    %waitforkey

end
%% Car 1 no-ett: Using mean as measurement

T = 0.1; % sample time 
% model 1 (turning)
f_point = @(st) [st(1)+T*st(3)*cos(st(4));...
                 st(2)+T*st(3)*sin(st(4));...
                 st(3);...
                 st(4) + T*st(5);...
                 st(5)];

   velCov1 = 0.5^2;  % velocity covariance
phiDotCov1 = 1^2; % turningrate covariance             
     subQ1 = diag([velCov1, phiDotCov1]);
gamma1 = [0 0 1 0 0;
          0 0 0 0 1]';
% MOTION COVARIANCE MATRIX 
Q1 = T*gamma1*subQ1*gamma1';    


x0 = [14.75 0.75 13 6.1 0]';
P0 = [ 0.0044    0.0002    0.0208   -0.0000   -0.0008;
       0.0002    0.0072    0.0054    0.0022   -0.0012;
       0.0208    0.0054    0.2942    0.0014   -0.0002;
      -0.0000    0.0022    0.0014    0.0012    0.0044;
      -0.0008   -0.0012   -0.0002    0.0044    0.1587];
       
rCov = 0.15;
h = @(st) [st(1); st(2)]; % Measurement function 

ukf = UKFsimple;
ukf.init(Q1, rCov, x0, P0);

for k = 1:length(data1)
    
    ukf.predictMoments(f_point);
    
    % Measurement 
    Mxy = mean(data1{k}(:,1:2));
    Z = Mxy';
    
    ukf.updateMoments({h},Z);
    
    x = ukf.upSt;
    P = ukf.upCov;
    
    car1no{k} = struct('st', x,'Cov',P);
    
    hold off
    plot(data1{k}(:,1), data1{k}(:,2),'x','Color', [.6 .6 .6]); axis equal; hold on    
    plot(Mxy(1), Mxy(2),'gs','MarkerSize',10)        
    plot(x(1), x(2), 'r*','MarkerSize',10)
    
    title(sprintf('%i/%i',[k length(data1)]))
    
    %waitforkey

end
 

%% Car 2 no-ett: Using mean as measurement

T = 0.1; % sample time 
% model 1 (turning)
f_point = @(st) [st(1)+T*st(3)*cos(st(4));...
                 st(2)+T*st(3)*sin(st(4));...
                 st(3);...
                 st(4) + T*st(5);...
                 st(5)];

   velCov1 = 0.5^2;  % velocity covariance
phiDotCov1 = 1^2; % turningrate covariance             
     subQ1 = diag([velCov1, phiDotCov1]);
gamma1 = [0 0 1 0 0;
          0 0 0 0 1]';
% MOTION COVARIANCE MATRIX 
Q1 = T*gamma1*subQ1*gamma1';    



x0 = [18.5 -13.4 4 pi/2 0]'; % Car 2

P0 = [ 0.0044    0.0002    0.0208   -0.0000   -0.0008;
       0.0002    0.0072    0.0054    0.0022   -0.0012;
       0.0208    0.0054    0.2942    0.0014   -0.0002;
      -0.0000    0.0022    0.0014    0.0012    0.0044;
      -0.0008   -0.0012   -0.0002    0.0044    0.1587];
       
rCov = 0.15;
h = @(st) [st(1); st(2)]; % Measurement function 

ukf = UKFsimple;
ukf.init(Q1, rCov, x0, P0);

data2 = car2;

for k = 1:length(data2)
    
    ukf.predictMoments(f_point);
    
    % Measurement 
    Mxy = mean(data2{k}(:,1:2));
    Z = Mxy';
    
    ukf.updateMoments({h},Z);
    
    x = ukf.upSt;
    P = ukf.upCov;
    
    car2no{k} = struct('st', x,'Cov',P);
    
    hold off
    plot(data2{k}(:,1), data2{k}(:,2),'x','Color', [.6 .6 .6]); axis equal; hold on    
    plot(Mxy(1), Mxy(2),'gs','MarkerSize',10)        
    plot(x(1), x(2), 'r*','MarkerSize',10)
    
    title(sprintf('%i/%i',[k length(data2)]))
    
    %waitforkey

end

%% Plot state and covariance for each car case, with ett and without 

data2 = car2;
fig = figure; fig.Position = [100 100 1600 800];
for k = 1:length(data2)
        
    hold off
    plot(data2{k}(:,1), data2{k}(:,2),'x','Color', [.6 .6 .6]); axis equal; hold on    
    
    
    st_ett = car2ett{k}.st;
    st_no  = car2no{k}.st;
    
    drawMyRide(st_ett,'b');
    
    plot(st_no(1), st_no(2), '*r','MarkerSize',10);
    
    xh = [st_no(1) st_no(1)+cos(st_no(4))];
    yh = [st_no(2) st_no(2)+sin(st_no(4))];
    plot(xh, yh,'-r')
    
    title(sprintf('%i/%i',[k length(data2)]))
    axis([14 30 -16 2])
    %waitforkey

end


%%



fig = figure; fig.Position = [100 100 1600 800];
for k = 1:length(data1)
        
    hold off
    plot(data1{k}(:,1), data1{k}(:,2),'x','Color', [.6 .6 .6]); axis equal; hold on    
    
    
    st_ett = car1ett{k}.st;
    st_no  = car1no{k}.st;
    
    drawMyRide(st_ett,'b');
    
    plot(st_no(1), st_no(2), '*r','MarkerSize',10);
    
    xh = [st_no(1) st_no(1)+cos(st_no(4))];
    yh = [st_no(2) st_no(2)+sin(st_no(4))];
    plot(xh, yh,'-r')
    
    title(sprintf('%i/%i',[k length(data1)]))
    %axis([14 30 -16 2])
    %waitforkey

end


%% Make comparison PLOTS! 

N1 = length(car1ett);
st1_ett = zeros(7,N1);
st1_no  = zeros(5,N1);

P1_ett  = zeros(7,7,N1);
P1_no   = zeros(5,5,N1);

for j = 1:N1
    st1_ett(:,j)  = car1ett{j}.st;
    P1_ett(:,:,j) = car1ett{j}.Cov;
    
    st1_no(:,j)  = car1no{j}.st;  
    P1_no(:,:,j) = car1no{j}.Cov;
end

N2 = length(car2ett);
st2_ett = zeros(7,N2);
st2_no  = zeros(5,N2);

P2_ett  = zeros(7,7,N2);
P2_no   = zeros(5,5,N2);

for j = 1:N2
    st2_ett(:,j)  = car2ett{j}.st;
    P2_ett(:,:,j) = car2ett{j}.Cov;
    
    st2_no(:,j)  = car2no{j}.st;    
    P2_no(:,:,j) = car2no{j}.Cov;

end



%%
state = 5;
fig = figure; fig.Position = [100 100 1600 800];
    subplot(2,1,1)
        
        p1 = plot(st1_ett(state,:),'--b'); hold on; axis tight
        p2 = plot(st1_no(state,:),'--r');        
        
        sigEtt = 3*sqrt(reshape(P1_ett(state,state,:), 1,N1));
        sigNo  = 3*sqrt(reshape(P1_no(state,state,:), 1,N1));
        
        shadeX = [1:N1, N1:-1:1];
        shadeY = [sigEtt+st1_ett(state,:) fliplr(-1.*sigEtt+st1_ett(state,:))];
        
        fl = fill(shadeX, shadeY, 'b');
        fl.FaceAlpha = 0.05;
        fl.EdgeColor = 'b';
        fl.EdgeAlpha = 0;
        
        shadeX = [1:N1, N1:-1:1];
        shadeY = [sigNo+st1_no(state,:) fliplr(-1.*sigNo+st1_no(state,:))];
        
        fl = fill(shadeX, shadeY, 'r');
        fl.FaceAlpha = 0.05;
        fl.EdgeColor = 'r';
        fl.EdgeAlpha = 0;
        
        plot(sigEtt+st1_ett(state,:),':b')
        plot(-1.*sigEtt+st1_ett(state,:),':b')
        plot(sigNo+st1_no(state,:),':r')
        plot(-1.*sigNo+st1_no(state,:),':r')
        
            t = title(sprintf('Car 1, state: %i',[state]));
            t.FontSize = 20;
            p1.LineWidth = 1;
            p2.LineWidth = 1;
    subplot(2,1,2)
        
        p1 = plot(st2_ett(state,:),'--b'); hold on; axis tight
        p2 = plot(st2_no(state,:),'--r');
        
        sigEtt = 3*sqrt(reshape(P2_ett(state,state,:), 1,N2));
        sigNo  = 3*sqrt(reshape(P2_no(state,state,:), 1,N2));
                
        shadeX = [1:N2, N2:-1:1];
        shadeY = [sigEtt+st2_ett(state,:) fliplr(-1.*sigEtt+st2_ett(state,:))];
        
        fl = fill(shadeX, shadeY, 'b');
        fl.FaceAlpha = 0.05;
        fl.EdgeColor = 'b';
        fl.EdgeAlpha = 0;
        
        shadeX = [1:N2, N2:-1:1];
        shadeY = [sigNo+st2_no(state,:) fliplr(-1.*sigNo+st2_no(state,:))];
        
        fl = fill(shadeX, shadeY, 'r');
        fl.FaceAlpha = 0.05;
        fl.EdgeColor = 'r';
        fl.EdgeAlpha = 0;
                
        plot(sigEtt+st2_ett(state,:),':b')
        plot(-1.*sigEtt+st2_ett(state,:),':b')
        plot(sigNo+st2_no(state,:),':r')
        plot(-1.*sigNo+st2_no(state,:),':r')
        
            t = title(sprintf('Car 2, state: %i',[state]));
            t.FontSize = 20;
            p1.LineWidth = 2;
            p2.LineWidth = 2;

%% Look at 

set(0,'defaulttextinterpreter','latex')
lim = [0.8, 1, 1.75, 1.25, 3];
fig = figure; fig.Position = [100 100 600 1000];

textLabels = {'$x_k$ $[m]$', '$y_k$ $[m]$', '$v_k$ $[m/s]$', '$\phi_k$  $[rad]$', '$\dot{\phi}_k$  $[rad/s]$'};

for k = 1:5
    state = k;
    subplot(5,1,k)
        hold on; axis auto; grid on; axis tight                        
        
        set(gca,'FontSize',12)
        
        sigEtt = 3*sqrt(reshape(P1_ett(state,state,:), 1, N1));
        sigNo  = 3*sqrt(reshape(P1_no(state,state,:), 1, N1));
        
        % The red one
        shadeX = [1:N1, N1:-1:1];
        shadeY = [sigNo fliplr(-1.*sigNo)];
        
        fl = fill(shadeX, shadeY, 'r');
        fl.FaceAlpha = 0.1;
        fl.EdgeColor = 'r';
        fl.EdgeAlpha = 0;
        
        % The blue one
        shadeX = [1:N1, N1:-1:1];
        shadeY = [sigEtt fliplr(-1.*sigEtt)];
        
        fl = fill(shadeX, shadeY, 'b');
        fl.FaceAlpha = 0.2;
        fl.EdgeColor = 'b';
        fl.EdgeAlpha = 0;        
        
        p1 = plot(sigEtt,':b');
        p2 = plot(-1.*sigEtt,':b');
        p3 = plot(sigNo,':r');
        p4 = plot(-1.*sigNo,':r');
        p1.LineWidth = 1.5;
        p2.LineWidth = 1.5;
        p3.LineWidth = 1.5;
        p4.LineWidth = 1.5;
        
        axis([1 N1 -lim(k) lim(k)])
        
        % Set ytick limits 
        myLim = linspace(-lim(k), lim(k), 5);
        set(gca,'ytick',myLim);
        
        % set ytick precision 
        lab = cell(1,5);
        for j = 1:5
            lab{1,j} = sprintf('%.3f',myLim(j));
        end        
        set(gca, 'YTickLabel', lab);                                
  

        yl = ylabel(textLabels{k},'fontsize',12);
        yl.FontSize = 20;        
        
        % adjust margins for subplot 
        pos = get(gca, 'Position');        
        pos(4) = 1/6.75;  %% This one here changes distance between vertical subplots
        pos(1) = 0.125;
        pos(3) = 0.825;
        set(gca, 'Position', pos)
        
        if k == 1
            t = title('Car 1 UKF $\pm3\sigma$-bounds');
            t.FontSize = 20;
        end
        
        if k == 5
            set(gca,'xtick',1:5:N1);
            xl = xlabel('\textbf{Time-step $k$}');
            xl.FontSize = 18;            
        else
            set(gca,'xtick',1:5:N1);
            set(gca,'xticklabel',[]);   
        end
        
end
        
        ax = gca;
        ax.XAxis.MinorTick = 'on';
        ax.XAxis.MinorTickValues = 1:N1

fig.PaperPositionMode = 'auto';
path = '~/Desktop/';
print(fig,strcat(path,sprintf('car1')),'-dpng','-r300','-opengl');
% 
% 
%% Calculate 
P1_ett
P1_no

P2_ett
P2_no
%%
rms = @(Xsquared) sqrt(sum(Xsquared)/length(Xsquared));

    for k = 1:5      
        alp = rms(reshape(P2_ett(k,k,:), 1, N2));
      fprintf('& %.3f ',[alp]);
    end






















