function y = gamma_2d_log(x)
%multivariate gamma function for the d=2 case

    y = (1/2)*log(pi) + gammaln(x) + gammaln(x-(1/2));
end