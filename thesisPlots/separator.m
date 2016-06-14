

rng(1337)
B = rand(200,2);

%plot(B(:,1), B(:,2),'rx'); hold on


f = @(x,y,m,k) [y-m-k*x > 0];

for j = 1:200
    if f(B(j,1),B(j,2),1,-1) &&  B(j,1) > 0.2 &&  B(j,2) > 0.2
        plot(B(j,1), B(j,2),'bo'); hold on
    else
        plot(B(j,1), B(j,2),'rx'); hold on
    end
end

axis([-.2 1.2 -.2 1.2])
%
x = 0:0.1:1;
y = -1.1*x + 1.1;


plot(x,y,'-k')

%% 


