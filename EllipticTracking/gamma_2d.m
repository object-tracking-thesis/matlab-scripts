function y = gamma_2d(x)
%multivariate gamma function for the d=2 case

    y = pi^(1/2) * gamma(x) * gamma(x - 1/2);
end