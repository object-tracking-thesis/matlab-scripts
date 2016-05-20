% A test to see how good the jacobians of mgpgen2 are by comparing them to
% numerical estimates

mgpgen2 = MGPgenerator2(1);

%%
%     x y v p pd w l
st = [1 2 3 4 5  6 7]';

p = 0.9;
K = 2;
j = 1:2;


jac1w = mgpgen2.mgpJac1_w;
fun1w = mgpgen2.mgpFunCorner1_w;

mgpgen2.evaluateJacobian(jac1w, st, K, 0, 0)

mgpgen2.evaluateFunction(fun1w, st, K, 0, 0);


d = 0.0001;
st1 = [1 2 3 4 5 6 7];

A0 = mgpgen2.evaluateFunction(fun1w, st, K, 0,0);
A1 = mgpgen2.evaluateFunction(fun1w, st1, K, 0,0);

deriv = (A1 - A0)./d;
deriv(2)
