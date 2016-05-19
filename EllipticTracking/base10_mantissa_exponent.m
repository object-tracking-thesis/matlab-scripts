function [mantissa, base10_exponent] = base10_mantissa_exponent(base,exponent)

act_exp = exponent*log10(abs(base));
base10_exponent = floor(act_exp);
mantissa = power(10,act_exp - base10_exponent);

if rem(exponent,2)==1 && sign(base)==-1
    mantissa = -mantissa;
end

return;
%credit to: https://stackoverflow.com/questions/25533237/how-to-compute-an-exponent-in-matlab-without-getting-inf