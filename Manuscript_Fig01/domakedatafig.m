% script to make and export figures for figure 1 of manuscript
clear; close all; clc;

% first simulation, on eupnea system
includec = 0;
total = 2040;
fig = zeros(1,11);

[param_e,out_e] = tabakrinzelcalcium('includec',includec,...
                                     'total',total,...
                                     'fig',fig);

% second simulation, just ca

% third simulation, both

% figure stuff
figure('Position',[5 5 1600 700])
spr = 3; % subplot rows
spc = 3; % subplot columns

% first row

% first column - traces
subplot(spr,spc,1)
plot(out_e.t,out_e.a)%,'b',out_e.t,out_e.s,'r',out_e.t,out_e.theta,'g')

% second column - nulclines
subplot(spr,spc,2)

thetaa  = param_e.thetaa;
ka      = param_e.ka;
lambdaa = param_e.lambdaa;
thetas  = param_e.thetas;
ks      = param_e.ks;
amax    = param_e.amax;
w       = param_e.w;
aa = 0:0.001:amax;
ss_lowtheta = zeros(length(aa),1);
ss_hightheta = zeros(length(aa),1);

% defining the two theta boundaries from the simulation
lowtheta  = min(out_e.theta);
hightheta = max(out_e.theta);
% loop through all a values and find the steady state s values
for i = 1:length(aa)
    a = aa(i);
    theta = lowtheta;
    fun = @(s) xinf(w*a*s-theta,thetaa,ka,lambdaa)-a;
    ss_lowtheta(i) = fzero(fun,0.5);
    theta = hightheta;
    fun = @(s) xinf(w*a*s-theta,thetaa,ka,lambdaa)-a;
    ss_hightheta(i) = fzero(fun,0.5);
end

plot(ss_lowtheta,aa,'r'); hold on
plot(ss_hightheta,aa,'r--');
plot(xinf(aa,thetas,ks,1),aa,'g');
plot(out_e.s,out_e.a,'b--')
axis([-0.1 1.2 0 amax+0.1])
legend("a'=0 w/ theta=2","a'=0 w/ theta=5","g'=0")

function y = xinf(x,theta,k,lambda)
y=lambda./(1+exp(4*(theta-x)/k));
end

