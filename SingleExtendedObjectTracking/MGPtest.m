% Define MGP functions for each car side

%    _<-3_
% ^ |     | ^
% | |  ? | |
% 2 |  x  | 4
%   |     |
%   |_____|
%    <-1
%
%       1    2    3    4      5         6    7
% st = [x_k, y_k, v_k, phi_k, phiDot_k, w_k, l_k]';

side1 = @(st) [...
               st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + 0.5*st(6)*cos(st(4) + pi/2);...
               st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + 0.5*st(6)*sin(st(4) + pi/2)
               ];

side2 = @(st) [...
               st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + st(6)*cos(st(4) + pi/2) + 0.5*st(7)*cos(st(4));...
               st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + st(6)*sin(st(4) + pi/2) + 0.5*st(7)*sin(st(4))
               ];

side3 = @(st) [...
               st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + 0.5*st(6)*cos(st(4) + pi/2) + st(7)*cos(st(4));...
               st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + 0.5*st(6)*sin(st(4) + pi/2) + st(7)*sin(st(4))
               ];
           
side4 = @(st) [...
               st(1) - 0.5*sqrt(st(6)^2 + st(7)^2)*cos(atan(st(6)/st(7)) + st(4)) + 0.5*st(7)*cos(st(4));...
               st(2) - 0.5*sqrt(st(6)^2 + st(7)^2)*sin(atan(st(6)/st(7)) + st(4)) + 0.5*st(7)*sin(st(4))
               ];
           
% Define box to do testing, this is in place of prediction
w = 1; l = 2;
Box = [0,0;...
       l,0;...
       l,w;...
       0,w;...
       0,0;...
       l/2,w/2];%x,y

offset = [0,0];
phi = pi/4;%2*pi*rand;
R = [cos(phi), -sin(phi); sin(phi), cos(phi)];

Box = Box*R' + repmat(offset,length(Box),1);

centerPoint = [l/2, w/2]*R' + offset;

st = [centerPoint(1), centerPoint(2), 0, phi, 0, w, l];

P1 = side1(st);
P2 = side2(st);
P3 = side3(st);
P4 = side4(st);

figure 
    grid on; hold on    
    plot(Box(:,1),Box(:,2),'-x') 
    plot(centerPoint(1),centerPoint(2),'ms')
    

    plot(P1(1),P1(2),'rx')
    plot(P2(1),P2(2),'rx')
    plot(P3(1),P3(2),'rx')
    plot(P4(1),P4(2),'rx')


    axis([-5 5 -5 5]); axis square








