figure;
for j=1:1
    pcshow(cleanedFrames{j})
    zoom(1.45)
    axis([-20 50 -40 20 -4 4])  
    pause(0.3)
end
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')
set(gca,'xtick',[-20:20:50]);         
set(gca,'ytick',[-40:20:20]);
set(gca,'ztick',[-4,0,4]);