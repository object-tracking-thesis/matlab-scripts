function color = class2color(clas)
% Takes a vector 'clas' containing numbers specifying class, and 
% generates these into RGB codes. 
%
% TODO: Currently pretty slow 
% 
    color = ones(length(clas),3);
    
    
    for k = 1:length(color)
        rng(clas(k))
        color(k,:) = uint8(rand(1,3)*255);
    end




end 