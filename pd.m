theta = 2;

ka = 2;
thetaa = 0;
w = 1;
lambdaa=10;
ks = -0.8;
thetas = 2.4;


a4s = linspace(0,10,1e4);
% a nulcline
s = (-ka*log(lambdaa./a4s -1)./4 + theta + thetaa)./(w*a4s);

s4a = linspace(0,1,1e4);
% s nulcine
a = (-ks*log(1./s4a-1))./4 + thetas;

figure
plot(s,a4s,'r'); hold on;
plot(s4a,a,'g');
legend("a'=0","s'=0")
axis([0 1 0 10])
