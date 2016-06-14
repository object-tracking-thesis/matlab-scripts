x = linspace(0,1,500);
y = -log(x);
subplot (1, 2, 1)
plotxy(x,y,1,'-log(x)')

x = linspace(0,1,500);
y = -log(1-x);
subplot (1, 2, 2)
plotxy(x,y,2,'-log(1-x)')