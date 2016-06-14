% % New tests for the PHD filter
% 
% % Test out CarTarget 
% 
% 
% carTarget = CarTarget;
% 
% load car1.mat
% load car2.mat
% 
% meas = car2adj;
% 
% N = 2;
% 
% plot(meas{N}(:,1), meas{N}(:,2), 'x','Color',[.6 .6 .6]); 
% axis equal; hold on
%     
%     %st0 = [14.75 0.75 13 6.1 0 1.8 4.5]'; % CAR 1
%      st0 = [18.5 -13.4 4 pi/2 0 1.6 3.5]';
%     %drawMyRide(st0,'b')
% %    
%     
% carTarget.init(st0, [], [], 1);
% 
% carTarget.predict();
% 
% predSt = carTarget.ukfmod.predSt;
%    drawMyRide(predSt, 'c')
%     
%     
% lik = carTarget.calcLikelihood(meas{N});
% yPred = carTarget.ukfmod.yPred;
% 
% for j = 1:length(yPred)
%    plot(yPred(1,j), yPred(2,j), 'mo') 
% end
% 
% drawMyRide(carTarget.ukfmod.upSt,'b')
% %%
% N = 3;
% plot(meas{N}(:,1), meas{N}(:,2), 'kx'); axis equal; hold on
% 
% carTarget.predict();
% 
% predSt = carTarget.ukfmod.predSt;
%     drawMyRide(predSt, 'c')
%     
%     
% lik = carTarget.calcLikelihood(meas{N});
% yPred = carTarget.ukfmod.yPred;
% 
% for j = 1:length(yPred)
%    plot(yPred(1,j), yPred(2,j), 'mo') 
% end
% 
% drawMyRide(carTarget.ukfmod.upSt,'b')

%% Let's test out the new PHD filter 

confirmedTargets = cell(1,186);

meas = load('~/Desktop/cm.mat');
meas = meas.cM;

% We define two targets, 
st1 = [14.75 0.75 13 6.1 0 1.8 4.5]'; % Car 1
st2 = [18.5 -13.4 4 pi/2 0 1.6 3.5]'; % Car 2

w = {1, 1};
brfs = {st1, st2};

phd = PHDinstance3(w, brfs);
fig = figure; hold on; axis equal
fig.Position = [50 50 1600 800];
%phd.K = 1;
for N = 1:186;

phd.predict();
hold off
% for j = 1:length(phd.componentStorage)
%     [st, cov] = phd.componentStorage(j).getState;
%     drawMyRide(st,'c'); hold on; axis equal
% end

% try
%     for m = 1:length(meas{N})
%         plot(meas{N}{m}(:,1), meas{N}{m}(:,2),'x'); hold on; axis equal
%     end
% catch me
% end

if ~isempty(meas{N}{1})
    phd.update(meas{N})
end

% for j = 1:length(phd.componentStorage)
%     [st, cov] = phd.componentStorage(j).getState;
%     drawMyRide(st,'b')
% end

T = phd.getBestRect(0.4);
tag = cell(1,length(T));
for j = 1:length(T)
   [st, cov]= T(j).getState();
   T(j).weight;
   drawMyRide(st,'r')
   tag{j} = T(j).getState();
end
confirmedTargets{1,N} = tag;
%title(sprintf('%i',N))
N;
if N == 10
    phd.K = 1;
end
% if N > 0
% waitforkey
% end
end


%% 

%mvnpdf()


plot(meas{38}{1}(:,1), meas{38}{1}(:,2),'kx'); hold on; axis equal
plot(meas{38}{2}(:,1), meas{38}{2}(:,2),'rx')









