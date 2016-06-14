%% Setup some parameters
input_layer_size  = 5;
hidden_layer_size = 10;
num_labels = 5;

%% Loading the data
fprintf('Loading Data ...\n')
% Load Training Data
%load lisa_data.mat
%load data/nn_clusters_kitti_static_crossing.mat
data = clusterObjects(randperm(size(clusterObjects,1)),:); %shuffling the data

X = data(:,1:input_layer_size);
y = data(:,input_layer_size+1);

m = size(X, 1);

[X mu_nn sigma_nn] = featureNormalize(X);
train = X(1:floor(0.6*m),:);
cv = X(ceil(0.6*m):floor(1.0*m),:);
test = X(ceil(0.8*m):m,:);
ytrain = y(1:floor(0.6*m),:);
ycv = y(ceil(0.6*m):floor(1.0*m),:);
ytest = y(ceil(0.8*m):m,:);

%% Initializing Parameters
initial_Theta1 = randInitializeWeights(input_layer_size, hidden_layer_size);
initial_Theta2 = randInitializeWeights(hidden_layer_size, num_labels);

% Unroll parameters
initial_nn_params = [initial_Theta1(:) ; initial_Theta2(:)];

%% Training the network
fprintf('\nTraining Neural Network... \n')

%parameters for the gradient descent
options = optimset('MaxIter', 200);
lambda = 1;

%shorthand function for the cost function
costFunction = @(p) nnCostFunction(p, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels, train, y, lambda);

%now, costFunction is a function that takes in only one argument (the neural network parameters)
[nn_params, cost] = fmincg(costFunction, initial_nn_params, options);

%obtain Theta1 and Theta2 back from nn_params
Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                 hidden_layer_size, (input_layer_size + 1));

Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end), ...
                 num_labels, (hidden_layer_size + 1));

fprintf('\npredicting %d clusters\n', size([cv; train; test],1));
tic
pred = predict(Theta1, Theta2, cv);
pred2 = predict(Theta1, Theta2, train);
pred3 = predict(Theta1, Theta2, test);
toc
%fprintf('\nCV Set Accuracy: %f\n', mean(double(pred == ycv)) * 100);
%fprintf('\nTrain Set Accuracy: %f\n', mean(double(pred2 == ytrain)) * 100);

%precision/recall
precs = [];
recas = [];
for i=2:num_labels
    PredY = [pred ycv]; %predictions and real outputs in one vector alongside each other
    tp = sum(PredY(:,1)==1 & PredY(:,2)==1);
    tn = sum(PredY(:,1)==i & PredY(:,2)==i);
    fp = sum(PredY(:,1)==1 & PredY(:,2)==i);
    fn = sum(PredY(:,1)==i & PredY(:,2)==1);
    prec = tp/(tp+fp);
    reca = tp/(tp+fn);
    precs = [precs prec];
    recas = [recas reca];
end

fprintf('\nCV Set Precision: %f\n', mean(precs));
fprintf('\nCV Set Recall: %f\n', mean(recas));
precs
recas
%PredY
