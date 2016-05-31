S = diag([1 1]);
v = [3 10 25];

figure;
for i = 1:9
    subplot(3,3,i)
    [x1,x2,x3] = threeSigmaOverGrid([0,0]',iwishrnd(S,v(ceil(i/3))));
    [xa1,xa2,xa3] = threeSigmaOverGrid([0,0]',S);
    plot(x3(1,:),x3(2,:))
    hold on
    plot(xa3(1,:),xa3(2,:),'--')
    hold on
    plot(0,0,'x')
    hold off
    legend('X','S')
    if (i > 6)
        xlabel('x')
    end
    if (mod(i,3) == 1)
        ylabel('y')
    end
end