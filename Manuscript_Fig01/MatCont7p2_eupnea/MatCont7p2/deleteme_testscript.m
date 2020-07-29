% quick a and s simulation w/ theta as parameter
close all;

theta=0.2;

t = 0:0.0001:60;
init = [1,0.2];

w=1;
lambdaa=10;
thetaa=0;
ka=2;
taua=0.15;
thetas=2.4;
ks=-0.8;
taus=0.75;

% % %

ainf = @(x) lambdaa./(1+exp(4* (thetaa-x)./ka));
sinf = @(a) 1./(1+exp(4*(thetas-a)./ks));

nsteps = length(t);
a = zeros(nsteps,1);
s = zeros(nsteps,1);
a(1) = init(1); s(1) = init(2);

dt = t(2)-t(1);

for i = 1:nsteps-1
    a(i+1) = a(i) + dt*(ainf(w*s(i).*a(i)-theta)-a(i))/taua;
    s(i+1) = s(i) + dt*(sinf(a(i))-s(i))/taus; 
end
    

plot(t,a)