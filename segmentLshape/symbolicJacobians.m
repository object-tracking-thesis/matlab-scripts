% Calculate partial derivaties for CASE 2

syms x y v phi phiDot W  L
% Function 1
fx1 = x - 0.5*sqrt(L^2 + W^2)*cos(atan(W/L) + phi);
fy1 = y - 0.5*sqrt(L^2 + W^2)*sin(atan(W/L) + phi);
 
fx1_dx = diff(fx1, x);
fx1_dy = diff(fx1, y);
fx1_dv = diff(fx1, v);
fx1_dphi = diff(fx1, phi);
fx1_dphiDot = diff(fx1, phiDot);
fx1_W = diff(fx1, W);
fx1_L = diff(fx1, L);

fy1_dx = diff(fy1, x);
fy1_dy = diff(fy1, y);
fy1_dv = diff(fy1, v);
fy1_dphi = diff(fy1, phi);
fy1_dphiDot = diff(fy1, phiDot);
fy1_W = diff(fy1, W);
fy1_L = diff(fy1, L);


J_F1 = [fx1_dx fx1_dy fx1_dv fx1_dphi fx1_dphiDot fx1_W fx1_L;...
        fy1_dx fy1_dy fy1_dv fy1_dphi fy1_dphiDot fy1_W fy1_L];
   
% Function 2
fx2 = fx1 + W*cos(phi + pi/2);
fy2 = fy1 + W*sin(phi + pi/2);

fx2_dx = diff(fx2, x);
fx2_dy = diff(fx2, y);
fx2_dv = diff(fx2, v);
fx2_dphi = diff(fx2, phi);
fx2_dphiDot = diff(fx2, phiDot);
fx2_W = diff(fx2, W);
fx2_L = diff(fx2, L);

fy2_dx = diff(fy2, x);
fy2_dy = diff(fy2, y);
fy2_dv = diff(fy2, v);
fy2_dphi = diff(fy2, phi);
fy2_dphiDot = diff(fy2, phiDot);
fy2_W = diff(fy2, W);
fy2_L = diff(fy2, L);

J_F2 = [fx2_dx fx2_dy fx2_dv fx2_dphi fx2_dphiDot fx2_W fx2_L;...
        fy2_dx fy2_dy fy2_dv fy2_dphi fy2_dphiDot fy2_W fy2_L];
     
% Function 3
fx3 = fx2 + L*cos(phi);
fy3 = fy2 + L*sin(phi);

fx3_dx = diff(fx3, x);
fx3_dy = diff(fx3, y);
fx3_dv = diff(fx3, v);
fx3_dphi = diff(fx3, phi);
fx3_dphiDot = diff(fx3, phiDot);
fx3_W = diff(fx3, W);
fx3_L = diff(fx3, L);

fy3_dx = diff(fy3, x);
fy3_dy = diff(fy3, y);
fy3_dv = diff(fy3, v);
fy3_dphi = diff(fy3, phi);
fy3_dphiDot = diff(fy3, phiDot);
fy3_W = diff(fy3, W);
fy3_L = diff(fy3, L);

J_F3 = [fx3_dx fx3_dy fx3_dv fx3_dphi fx3_dphiDot fx3_W fx3_L;...
        fy3_dx fy3_dy fy3_dv fy3_dphi fy3_dphiDot fy3_W fy3_L];
    
% tic    
% x = 1; y = 1; v = 1; phi = 1; phiDot = 1; W = 1; L = 1;
% eval(J_F1)
% eval(J_F2)
% eval(J_F3)
% toc   
   
%% Check how well the functions work 

x = 5; y = 5; phi = pi/4; % 45 deg angle 
W = 1; L = 2;

f = figure; grid on; hold on; xlabel('x'); ylabel('y');
    plot(x,y,'ko')
    plot([x x]+[0 cos(phi)], [y y]+[0 sin(phi)],'k-')
    axis([0 10 0 10]); axis square

    x1 = eval(fx1);
    y1 = eval(fy1);
    x2 = eval(fx2);
    y2 = eval(fy2);
    x3 = eval(fx3);
    y3 = eval(fy3);
    
    plot(x1,y1,'or')
    plot(x2,y2,'ob')
    plot(x3,y3,'og')

    plot([x1 x2 x3],[y1 y2 y3],'-c')








   
   
   
   