% Gating function for entire car. p is ratio of how much the car should be
% extended, i.e. p = 1.2 will result in 20% larger sides 
%
% bool = carGating(cluster, predictedState, p)
%
%


function bool = carGating(cluster, predictedState, p)

%state = [0 1 0 pi/6 0 1.8 4.5]';
    
% Extend state with p

%p = 1.5;

    st2 = predictedState;
    
    st2(end) = st2(end)*p;
    st2(end-1) = st2(end-1)*p;
    
    % Get lower point in rectangle (corner i)
    A = st2(1:2) - 0.5*sqrt(st2(6)^2 +  st2(7)^2)*[cos(atan(st2(6)/st2(7))+st2(4));...
                                                   sin(atan(st2(6)/st2(7))+st2(4))];                                                    
    A = A';                                                    
    
    AB = st2(7)*[cos(st2(4)) sin(st2(4))];

    AD = st2(6)*[cos(st2(4) + pi/2) sin(st2(4) + pi/2)];
    
    try 
        Mz = mean(cluster(:,1:2));
        bool = isInRectangle(Mz, A, AB, AD);
    catch me
        bool = 0;        
    end            
end
    
    
function bool = isInRectangle(Mz, A, AB, AD)
    % For rectangle 1 (upper)
    AM = Mz - A;

    bool = ((0 < dot(AM,AB)) && (dot(AM,AB)< dot(AB,AB))) && ((0 < dot(AM,AD)) && (dot(AM,AD)< dot(AD,AD)));
end
    
    
    
    
    
    
    
    
    
    
    
    
    