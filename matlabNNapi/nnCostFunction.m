function [J grad] = nnCostFunction(nn_params, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels, ...
                                   X, y, lambda)
%NNCOSTFUNCTION Implements the neural network cost function for a two layer
%neural network which performs classification
%   [J grad] = NNCOSTFUNCTON(nn_params, hidden_layer_size, num_labels, ...
%   X, y, lambda) computes the cost and gradient of the neural network. The
%   parameters for the neural network are "unrolled" into the vector
%   nn_params and need to be converted back into the weight matrices. 
% 
%   The returned parameter grad should be a "unrolled" vector of the
%   partial derivatives of the neural network.
%

% Reshape nn_params back into the parameters Theta1 and Theta2, the weight matrices
% for our 2 layer neural network
Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                 hidden_layer_size, (input_layer_size + 1));

Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end), ...
                 num_labels, (hidden_layer_size + 1));

% Setup some useful variables
m = size(X, 1);
         
J = 0;
Theta1_grad = zeros(size(Theta1));
Theta2_grad = zeros(size(Theta2));

%reforming y from digits to a classification vector Y
Y = zeros(m,num_labels);
for i = 1:m
	Y(i,y(i)) = 1;
end

%adding bias unit to X
X = [ones(m,1), X];

%feedforward from input layer to hidden layer
a2 = sigmoid(Theta1*X')';
%adding bias unit to a2
a2 = [ones(size(a2,1),1), a2];

%feedforward from hidden layer to output layer
h = sigmoid(Theta2*a2')';

%computing the total cost without regularization
for i = 1:size(h,1)
	J = J + (1/m)*sum(-Y(i,:).*log(h(i,:))-((ones(1,num_labels)-Y(i,:)).*log(ones(1,num_labels)-h(i,:))));
end

%adding the regularization term
%don't regularize the bias unit theta terms
theta1_1 = Theta1(:,2:end);
theta2_1 = Theta2(:,2:end);
J = J + (lambda/(2*m))*(sum(theta1_1(:).^2)+sum(theta2_1(:).^2));


%computing the gradient terms
delta1 = zeros(size(Theta1));
delta2 = zeros(size(Theta2));

for i = 1:m
	%adding bias unit to X
	xb = X(i,:);

	%feedforward from input layer to hidden layer
	z2 = Theta1*xb';
	a2 = sigmoid(z2)';
	%adding bias unit to a2
	a2 = [1, a2];

	%feedforward from hidden layer to output layer
	h = sigmoid(Theta2*a2')';

	%backprop
	%error for output layer
	d3 = h-Y(i,:);

	%error for hidden layer
	d2 = (Theta2'*d3').*sigmoidGradient([1; z2]);

	%accumulate the gradient
	delta1 = delta1 + d2(2:end)*xb;
 	delta2 = delta2 + d3(:)*a2;	
end

D1 = (1/m).*delta1;
D2 = (1/m).*delta2;

%Regularization of the gradients
for j = 2:size(Theta1,2)
	D1(:,j) = D1(:,j) + (lambda/m)*Theta1(:,j);
end
for j = 2:size(Theta2,2)
	D2(:,j) = D2(:,j) + (lambda/m)*Theta2(:,j);
end

Theta1_grad = D1;
Theta2_grad = D2;

% Unroll gradients
grad = [Theta1_grad(:) ; Theta2_grad(:)];

end
