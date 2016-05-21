% Draw a steering wheel 

for k = 1:length(stateStorage)
hold off
r = 1:0.01:1.1;
phi = 0:0.01:2*pi;

for j = r
   plot(j*cos(phi), j*sin(phi),'b'); hold on   
end
axis equal

phiDot = stateStorage{k}(5);
l = 0:0.1:1.05;
plot(cos(phiDot + pi/2).*l, sin(phiDot + pi/2).*l,'-r','MarkerFaceColor','r','LineWidth',16)

pause(0.1)

end