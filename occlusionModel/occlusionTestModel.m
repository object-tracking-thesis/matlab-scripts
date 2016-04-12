% Test to write occlusion model

% Generate random 2D point clouds
N = 15;

Q  = 0.01.*[1 0.2;0.2 1];
a = -3;
b = 3;


K = 10;
PC = cell(1,K);
rng(1335)

for k = 1:K
    N = randi(100);
    r = a*ones(2,1) + (b-a).*rand(2,1);
    PC{k} = mvnrnd(r, Q,N);
end
close all

pHandle = zeros(1,K);
strings = cell(1,K);
Pd = ones(K,1);
%%
figure('Position',[600 100 1000 800])

PC{1} = PC{1} + repmat(mvnrnd([-1 -2],Q,1),length(PC{1}),1);
for H = 1:1
    clf
    for k = 1:K
        PC{k} = PC{k} + repmat(mvnrnd([0 0],Q,1),length(PC{k}),1);
    end

    Pd = getDetectionProbability(PC);
    hold on
    p1 = plot(0,0,'ko','MarkerFaceColor','k');
    legend(p1,'EGO Vehicle')
    for k = 1:K
        
        plot(PC{k}(:,1),PC{k}(:,2),'.','MarkerSize',20)
        idx = getOcclusionPoints2D(PC{k});
        plot([0 PC{k}(idx(1),1)],[0 PC{k}(idx(1),2)],'k')
        plot([0 PC{k}(idx(2),1)],[0 PC{k}(idx(2),2)],'k')
        plot(PC{k}(:,1),PC{k}(:,2),'.','MarkerSize',20,'MarkerEdgeColor',[1 1 1]-repmat(Pd(k),1,3));
        strings{k} = sprintf('PointCloud %d Pd: %.2f',[k Pd(k)]);
        m = mean(PC{k});
        text(m(1)*1.1,m(2)*1.1,strings{k})
    end
    axis([-5 5 -5 5])
    pause(1)
    %waitforbuttonpress
end


%% Test out funcitons

x = 0:0.01:30;
a = 0;
b = 15;
slope = 0.8;
lowerLimit = 0.05;

y1 = lowerLimit/2 + 0.9./exp(slope.*(x-a));
y2 = lowerLimit/2 + 0.9./exp(-slope.*(x-b));

%plot(x,y1,x,y2)
hold on

plot(x,y1+y2,'k')

axis([0 20 0 1])

%%

a = 0;
b = 10;

I = mat2gray([a b])

alpha = 2;
L = 0.1;
x = 0:0.1:10;
f = @(x) L + (1-L)*(exp(-alpha*(x-a)) + exp(alpha*(x-b)));
y = f(x)

plot(x,y)
axis([0 10 0 1])




