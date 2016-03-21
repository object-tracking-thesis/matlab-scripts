function plotCylinder(centerPoint, radius, groundPlaneZ, height)

 phi = linspace(0, 2*pi, 10);
 xCirc = radius*cos(phi) + centerPoint(1);
 yCirc = radius*sin(phi) + centerPoint(2);
 
 X = [xCirc;xCirc];
 Y = [yCirc;yCirc];
 Z = [groundPlaneZ*ones(size(xCirc));...
     (groundPlaneZ+height)*ones(size(xCirc))];
                        
surf(X,Y,Z);
end