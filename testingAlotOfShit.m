

meas = load('~/Desktop/cm.mat');
meas = meas.cM;


% for k = 1:length(meas)
%     hold off
%    for j = 1:length(meas{k})
%       A = meas{k}{j};
%       hold on
%       plot(A(:,1), A(:,2), 'kx')
%    end
%    axis equal
%    waitforkey
% 
% end

%%




weights = {1, 1};

x1 = [14.5 0.75 4 6 0 1.8 4.7]';

%x2 = [18.4 -13.5 4 pi/2 0 1.8 4.7]';
x2 = [18.3, -14.5, 4.4, 1.51, 0, 1.8, 4]';       

%x3 = [24, -15, 4.4, 0, 0, 1.8, 2]';       


means = {x1, x2};
N = 22;
%     plot(meas{N}{1}(:,1), meas{N}{1}(:,2), 'kx'); hold on
% 
%     drawMyRide(x1)
%     drawMyRide(x2); axis equal


%%

phd = PHDinstance2(weights, means);
%     phd.predict();
%     
%     
%     
%         plot(meas{1}{1}(:,1), meas{1}{1}(:,2), 'kx'); hold on
% 
%     drawMyRide(phd.componentStorage(1).imm.getState)
%     
%     axis equal
%     
%%
fig = figure;    
for k = 1:186
    
    hold off
    plot3(0,0,0 ,'sg')
    hold on
    try
        plot3(meas{k}{1}(:,1), meas{k}{1}(:,2), meas{k}{1}(:,3), 'kx')
    catch me
    end    
    phd.predict();
    
    if ~isempty(meas{k}{1})
        phd.update(meas{k});
    end
    
    if k > 50
        phd.K = 1;
    end
    
    if k > 80
         phd.K = 100;
    end

    
    try
        plot3(meas{k}{2}(:,1), meas{k}{2}(:,2), meas{k}{2}(:,3), 'bx')
    catch me
    end
    
     A = phd.getBestRect(0.999);
     st = cell(1,1);
     for j = 1:length(A)
         A(j).weight
         [st, cov, w1, w2] = A(j).getState;          
         drawMyRideCube(st, -1.8,0, 'r')
         text(double(st(1)), double(st(2))+1, 0, sprintf('w1: %.2f, w2: %.2f',[w1 w2]));
     end    
    title(sprintf('k=%i',k))
    axis([0 42 -16 2 -2 2])
    axis equal
    pause(0.1)
    
end


%%




