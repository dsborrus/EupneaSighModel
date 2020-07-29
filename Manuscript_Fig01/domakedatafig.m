%% Manuscript Figure #01 - Introduction to the model
% script to make and export figures for figure 1 of manuscript
clear; close all; clc;

% figure stuff
figure('Position',[5 5 1600 700])
%figure(9)
spr = 3; % subplot rows
spc = 3; % subplot columns

%% first row, simulation, on eupnea system
includec = 0;
total = 2100;
fig = zeros(1,11); fig(2)=0;

[param_e,out_e] = tabakrinzelcalcium('includec',includec,...
                                     'total',total,...
                                     'fig',fig);
                                 
% first column - traces
subplot(spr,spc,1)
plot(out_e.t,out_e.a)%,'b',out_e.t,out_e.s,'r',out_e.t,out_e.theta,'g')
xlabel('t'); ylabel('a')

% second column - nulclines
subplot(spr,spc,2)

thetaa  = param_e.thetaa;
ka      = param_e.ka;
lambdaa = param_e.lambdaa;
thetas  = param_e.thetas;
ks      = param_e.ks;
amax    = param_e.amax;
w       = param_e.w;

% defining the two theta boundaries from the simulation
lowtheta  = min(out_e.theta);
hightheta = max(out_e.theta);

a4s = linspace(0,lambdaa,1e4);
% a nulcline
s_lowtheta  = (-ka*log(lambdaa./a4s -1)./4 + lowtheta + thetaa)./(w*a4s);
s_hightheta = (-ka*log(lambdaa./a4s -1)./4 + hightheta + thetaa)./(w*a4s);
s4a = linspace(0,1,1e4);
% s nulcine
a = (-ks*log(1./s4a-1))./4 + thetas;

plot(s_lowtheta,a4s,'r'); hold on
plot(s_hightheta,a4s,'r--');
plot(xinf(a4s,thetas,ks,1),a4s,'g');
plot(out_e.s,out_e.a,'b--');
axis([-0.1 1.2 0 amax+0.1])
legend("a'=0 w/ theta=2","a'=0 w/ theta=5","s'=0")
xlabel('s'); ylabel('a');
%comet(out_e.s,out_e.a)%,'b--')    

% Third column - theta bifurcation

subplot(spr,spc,3); hold on;

%% Second row - just ca system

includec = 1;
includes = 0;
includetheta = 0;
total = 2100;

[param_s,out_s] = tabakrinzelcalcium('includec',includec,...
                                     'includes',includes,...
                                     'includetheta',includetheta,...
                                     'total',total,...
                                     'fig',fig);
                                 
% first col - simple simulation
subplot(spr,spc,4); hold on;
plot(out_s.t,out_s.c,'k')
xlabel('t'); ylabel('c')

% second col - phase plane
subplot(spr,spc,5); hold on;

v1=20;
v2=0.25;
v3=60;
k3=0.3;
n3=2;
lambda=0.15;
thetam=0.25;
km=0.04;
thetah=0.3;
kh=-0.06;
j0=0.009;
v4=0.4;
k4=0.3;
n4=4;

j0 = param_s.jin0;
v4 = param_s.v4;
k4 = param_s.k4;

cc = 0:0.001:0.4;

for i=1:length(cc)
    cnull(i)=fzero(@(x) (v2+v1*finf(cc(i),thetam,km,thetah,kh))*((x-cc(i))/lambda-cc(i))-v3*cc(i)^n3/(k3^n3+cc(i)^n3)+j0-((v4*cc(i)^4)/(k4+cc(i)^4)),0.3);
end

% c nulc line
plot(cc,cnull,'m','linewidth',2);
% ct nulc line
plot([((j0*k4^4)/(v4-j0))^(1/4) ((j0*k4^4)/(v4-j0))^(1/4)], [0 max(cnull)],'k','linewidth',2 )
% trajectory
plot(out_s.c,out_s.ct,'b')
xlabel('c'); ylabel('ct')

% third column - bifurcation
subplot(spr,spc,6); hold on;

eq = load('../Fig03CalciumClosedCellModelVer02/MatCont7p2/Systems/sigh_bif_control/diagram/Step2_findequill_hopfs.mat');
eq_filler = load('../Fig03CalciumClosedCellModelVer02/MatCont7p2/Systems/sigh_bif_control/diagram/EP_EP(2).mat');
lc = load('../Fig03CalciumClosedCellModelVer02/MatCont7p2/Systems/sigh_bif_control/diagram/Step3_LimitCycle.mat');
    
% equill line
j0_eq = eq.x(3,:);
css = eq.x(1,:);

hopfs = [16,39];

% equll filler line
j0_eqfill = eq_filler.x(3,:);
css_eqfill = eq_filler.x(1,:);

% limit cycle
c_lc = lc.x(1:2:end-2,:);
j0_lc = lc.x(324,:);
period = lc.x(323,:);


% figures
plot(j0_eq(1:hopfs(1)),css(1:hopfs(1)),'k'); hold on;
plot(j0_eq(hopfs(1):hopfs(2)),css(hopfs(1):hopfs(2)),'k--');
plot(j0_eq(hopfs(2):end),css(hopfs(2):end),'k');
plot(j0_eq(hopfs),css(hopfs),'bo')
plot(j0_eqfill,css_eqfill,'k')
xlim([0 0.2])
% idk
temp = 475;
plot(j0_lc(1:temp),max(c_lc(:,1:temp)),'r');
plot(j0_lc(1:temp),min(c_lc(:,1:temp)),'r');
xlabel('j0')
ylabel('c ss')
% 
% subplot(2,1,2)
% plot(j0_lc(150:temp-42),period(150:temp-42),'linewidth',2)
% xlim([0 0.2])
% xlabel('jin0')
% ylabel('period')

% third simulation, both

function y = xinf(x,theta,k,lambda)
y=lambda./(1+exp(4*(theta-x)/k));
end

function f = finf(x,theta1,k1,theta2,k2)
    f=1./(1+exp((theta1-x)/k1))./(1+exp((theta2-x)/k2));
end
