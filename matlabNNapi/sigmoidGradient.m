function g = sigmoidGradient(z)
%SIGMOIDGRADIENT returns the gradient of the sigmoid function
%evaluated at z
g = zeros(size(z));
g = sigmoid(z).*(ones(size(z)) - sigmoid(z));

end
