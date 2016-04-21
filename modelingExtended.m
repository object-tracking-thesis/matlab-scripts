% Look at modelling the density of an extended target 

L = 4.58; % meters
B = 1.80; % meters
 hold on
s1 = [0 0];
P1 = 0.01*[0.1 0;
      0 1];
N = 200;
[X1,X2] = meshgrid(linspace(-2,8,N)', linspace(-2,8,N)');
X = [X1(:) X2(:)];
p = mvnpdf(X, s1, P1);
sups = surf(X1,X2,reshape(p,N,N));
sups.FaceAlpha = 0.2;
sups.EdgeAlpha = 0.2;

s1 = [0 L];
P1 = 0.01*[0.1 0;
      0 1];
[X1,X2] = meshgrid(linspace(-2,8,N)', linspace(-2,8,N)');
X = [X1(:) X2(:)];
p = mvnpdf(X, s1, P1);
sups = surf(X1,X2,reshape(p,N,N));
sups.FaceAlpha = 0.2;
sups.EdgeAlpha = 0.2;

X = [0 0 B B 0;
     0 L L 0 0]';
plot(X(:,1),X(:,2),'r','LineWidth',2)
axis([-2 8 -2 8])


%% 
N = 1500;
Z = mvnrnd(ones(N,1),ones(N,N),1)';
tic 
prob = mvnpdf(Z,ones(N,1),diag([ones(1,N)]));
toc
%%
prod = 1;
tic 
for j = 1:N
    prod = prod * mvnpdf(Z(1),1,1);
end
toc