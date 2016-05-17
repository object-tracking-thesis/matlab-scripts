function giw_comp = giwComp(mu,P,v,S,weight,index)
%a gaussian inverse-wishart component
%parametrized by:
% mu    - mean of the gaussian
% P     - covariance of the gaussian
% v     - scalar degrees of freedom for the IW
% S     - DxD positive definite matrix for the IW
% weight- weight of the component
% index - index of the component

    giw_comp = struct('mu',mu,'P',P,'v',v,'S',S,'weight',weight,'index',index);
end