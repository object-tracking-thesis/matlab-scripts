%% Generate multiple tracks, 2D xy
%rng(2) % 20,36 is good 
T = 1;
A = [1 T 0 0;
    0 1 0 0;
    0 0 1 T;
    0 0 0 1];

mu_s = [0 0 0 0];

Q_s = 0.05.*[0 0 0 0;
             0 1 0 0;
             0 0 0 0;
             0 0 0 1];

H = [1 0 0 0;
     0 0 1 0];

mu_m = [0 0];

Q_m = 0.025.*[1 0;
              0 1];

x1_init = [1 2 -10 0]';
x2_init = [20 0 10 -2]';
x3_init = [5 1 -10 1]';
x4_init = [14 0 -18 2]';
% Generate State & Measurement Sequence for two targets

N = 100;
vLimit = round(N/2);
% x1 & x2 always present. x3 is birthed halfway through, x4 dies halfway through 
x1_state = ones(4,N);
x2_state = ones(4,N);

x3_state = ones(4,vLimit);
x4_state = ones(4,vLimit);

x1_m = ones(3,N); % Bottom row stands for detected or not 
x2_m = ones(3,N); % Bottom row stands for detected or not 
x3_m = ones(3,vLimit); % Bottom row stands for detected or not 
x4_m = ones(3,vLimit); % Bottom row stands for detected or not 


% Set probability of detection

Pd = 0.75;
% Set states 
for j = 1:N
    if j == 1
        x1_state(:,j) = A*x1_init + mvnrnd(mu_s,Q_s)';
        x2_state(:,j) = A*x2_init + mvnrnd(mu_s,Q_s)';
        
        x1_m(1:2,j) = H*x1_state(:,j) + mvnrnd(mu_m,Q_m)';
        x2_m(1:2,j) = H*x2_state(:,j) + mvnrnd(mu_m,Q_m)';
                
        x1_m(3,j) = binornd(1,Pd);
        x2_m(3,j) = binornd(1,Pd);
    else
        x1_state(:,j) = A*x1_state(:,j-1) + mvnrnd(mu_s,Q_s)';
        x2_state(:,j) = A*x2_state(:,j-1) + mvnrnd(mu_s,Q_s)';
        
        x1_m(1:2,j) = H*x1_state(:,j) + mvnrnd(mu_m,Q_m)';
        x2_m(1:2,j) = H*x2_state(:,j) + mvnrnd(mu_m,Q_m)';
                
        x1_m(3,j) = binornd(1,Pd);
        x2_m(3,j) = binornd(1,Pd);
    end
end

% States that die or get birthed 
vLimit = round(N/2);
for j = 1:vLimit
    if j == 1
        x3_state(:,j) = A*x3_init + mvnrnd(mu_s,Q_s)';
        x4_state(:,j) = A*x4_init + mvnrnd(mu_s,Q_s)';
        
        x3_m(1:2,j) = H*x3_state(:,j) + mvnrnd(mu_m,Q_m)';
        x4_m(1:2,j) = H*x4_state(:,j) + mvnrnd(mu_m,Q_m)';
                
        x3_m(3,j) = binornd(1,Pd);
        x4_m(3,j) = binornd(1,Pd);
    else
        x3_state(:,j) = A*x3_state(:,j-1) + mvnrnd(mu_s,Q_s)';
        x4_state(:,j) = A*x4_state(:,j-1) + mvnrnd(mu_s,Q_s)';
        
        x3_m(1:2,j) = H*x3_state(:,j) + mvnrnd(mu_m,Q_m)';
        x4_m(1:2,j) = H*x4_state(:,j) + mvnrnd(mu_m,Q_m)';
                
        x3_m(3,j) = binornd(1,Pd);
        x4_m(3,j) = binornd(1,Pd);
    end
end

%% Generate Clutter

% Set density rate 
be = 5e-2;
% Set Volume/Area 
Vl = 36*42;

% Draw poission value, use this to generate uniformly distributed xy values

clutterMeas = cell(1,N);
xa = -1;
xb = 45;
ya = -40;
yb = 12;
for k = 1:N
   poissValue = poissrnd(be*Vl);
   Xc = xa + (xb - xa)*rand(1,poissValue);
   Yc = ya + (yb - ya)*rand(1,poissValue);
   clutterMeas{1,k} = [Xc;Yc];
end



%% (New) Plot the trajectories

X1 = [x1_init x1_state];
X2 = [x2_init x2_state];
X3 = [x3_init x3_state];
X4 = [x4_init x4_state];

X1m = [[0 0 0]' x1_m];
X2m = [[0 0 0]' x2_m];
X3m = [[0 0 0]' x3_m];
X4m = [[0 0 0]' x4_m];

if 0
    button = 'Yes';
    h1 = figure('Position',[100, 400, 1920*0.7 1080*0.7]);
    
    k = 1;
    checker = 0;
    vLimit = vLimit +1;
    while h1.isvalid == 1 && strcmp('Yes',button)
        %for k = 1:N+1
        while k > 0 && k < N +2 && h1.isvalid
            hold on
            
            if k > 1
                % Plot Clutter
                hold off
                plot(clutterMeas{k-1}(1,:),clutterMeas{k-1}(2,:),'x','Color',[0.3 0.3 0.3])
                hold on
                %plot(X1m(1,2:k),X1m(2,2:k),'xb')
                plotMeasurement(X1m(1,2:k),X1m(2,2:k),X1m(3,k),'b')
                plotMeasurement(X2m(1,2:k),X2m(2,2:k),X2m(3,k),'r')
            end
            plot(X1(1,1:k),X1(3,1:k),'ob')
            plot(X2(1,1:k),X2(3,1:k),'or')
            
            if k <= vLimit
                plot(X3(1,1:k),X3(3,1:k),'og')
                plotMeasurement(X3m(1,2:k),X3m(2,2:k),X3m(3,k),'g')
                plotEllip(X3(1:2:3,k), Q_m, [1 2])
            else
                a = k-vLimit;
                plot(X4(1,1:a),X4(3,1:a),'om')
                plot(X4m(1,2:a),X4m(2,2:a),'xm')
                plotEllip(X4(1:2:3,a), Q_m, [1 2])
            end
            
            axis([-1 45 -40 12])
            xlabel('X')
            ylabel('Y')
            title(sprintf('TimeStep: %d/%d',[k-1, N]))
            
            plotEllip(X1(1:2:3,k), Q_m, [1 2])
            plotEllip(X2(1:2:3,k), Q_m, [1 2])
            
            if k < N+1
                a = 1+k;
                [xf, yf]=ds2nfu(X1(1,k:a),X1(3,k:a));
                xf(xf < 0) = 0;
                yf(yf < 0) = 0;
                annotation('line', xf,yf);
                
                [xf, yf]=ds2nfu(X2(1,k:a),X2(3,k:a));
                xf(xf < 0) = 0;
                yf(yf < 0) = 0;
                annotation('line', xf,yf)
                
                if k <= vLimit-1
                    [xf, yf]=ds2nfu(X3(1,k:a),X3(3,k:a));
                    xf(xf < 0) = 0;
                    yf(yf < 0) = 0;
                    annotation('line', xf,yf)
                elseif k > vLimit
                    [xf, yf]=ds2nfu(X4(1,k-vLimit:a-vLimit),X4(3,k-vLimit:a-vLimit));
                    xf(xf < 0) = 0;
                    yf(yf < 0) = 0;
                    annotation('line', xf,yf)
                end
            end
            
            if checker == 10
                pause(0.750)
            end
            
            %checker
            if  checker ~= 10
                checker = mydialog;
                if checker ~= 10
                    k = k + checker;
                end
            else
                k = k +1;
            end
            
        end
        
        if h1.isvalid
            button = questdlg('Replay Sequence?','MHT Player', 'Yes','No','Yes');
            if strcmp(button,'Yes')
                clf;
                k = 1;
                checker = 0;
            else
                close(h1);
            end
        end
    end
end





%% Old plot
st = 1;
sp = N;
h1 = figure('Position',[50 50 1920*0.8 1080*0.8]);
    hold on
    plot(x1_state(1,st:sp),x1_state(3,st:sp),'ob')
    plot(x1_init(1),x1_init(3),'ob','MarkerFaceColor','blue')
    plot(x1_m(1,st:sp),x1_m(2,st:sp),'xb')

    plot(x2_state(1,st:sp),x2_state(3,st:sp),'or')
    plot(x2_init(1),x2_init(3),'or','MarkerFaceColor','red')
    plot(x2_m(1,st:sp),x2_m(2,st:sp),'xr')

    xlabel('X')
    ylabel('Y')

    X1 = [x1_init x1_state];
    X2 = [x2_init x2_state];
%     for k = 1:length(X1)-1
%         a = 1+k;
%         [xf, yf]=ds2nfu(X1(1,k:a),X1(3,k:a));
%         annotation('arrow', xf,yf)
%         [xf, yf]=ds2nfu(X2(1,k:a),X2(3,k:a));
%         annotation('arrow', xf,yf)
%     end




