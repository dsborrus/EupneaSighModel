function [ param, out, mmc ] = calciumsystem(varargin)

% system parameters
param.trans = 2000;
param.total = 2200;
param.dt=0.001;
param.seed=-1;
param.ct_iei_binwidth=5;

% model parameters 
param.ct_thresh=1.2; % uM; threshold crossing for total calcium
param.ct_ieimin=0;

param.v1=20;
param.v2=0.25;
param.v3=60;
param.k3=0.3;
param.n3=2;
param.lambda=0.15;
param.thetam=0.25;
param.km=0.04;
param.thetah=0.3;
param.kh=-0.06;

param.jin0=0.009; % 0.02 0.018 0
param.jin1=0; % 0.02 0.04 0.2
param.v4=0.4;
param.k4=0.3;
param.n4=4;

param.thetac=0.35; % 0.35
param.kc=0.05; % 0.05

param.init_c=0.1;
param.init_ct=2;

param.trjwnd = 200;

param.compactout = 0; % in order to not return actual simulation results
                      % i.e. (t,a,s,theta,c,ct) so save on space during
                      % pstudy
 
if nargin>0
    disp('Non-default parameters: ')
    for i=1:nargin/2
        param.(varargin{2*i-1})=varargin{2*i};
        disp(['   ' varargin{2*i-1} ' = ' num2str(varargin{2*i}) ])
    end
end

trans=param.trans;
total=param.total;
dt=param.dt;
seed=param.seed;

ct_iei_binwidth=param.ct_iei_binwidth;

ct_thresh=param.ct_thresh;
ct_ieimin=param.ct_ieimin;

v1=param.v1;
v2=param.v2;
v3=param.v3;
k3=param.k3;
n3=param.n3;
lambda=param.lambda;
thetam=param.thetam;
km=param.km;
thetah=param.thetah;
kh=param.kh;

jin0=param.jin0;
jin1=param.jin1;
v4=param.v4;
k4=param.k4;
n4=param.n4;

thetac=param.thetac;
kc=param.kc;

init_c=param.init_c;
init_ct=param.init_ct;

close all
if seed==-1, rng('shuffle'); else rng(seed); end

%% integration of ODEs
t = [0:dt:total]; % for integration of odes; times includes transient
[c, ct] = deal(zeros(size(t)));
c(1)=init_c; ct(1)=init_ct;
ct_cross = [];

for i=2:length(t)

    a=0;
    jpm = jin0 + jin1*a - v4*c(i-1)^n4/(k4^n4+c(i-1)^n4);
    c(i)=c(i-1)+dt*(  ( (v2+v1*finf(c(i-1),thetam,km,thetah,kh))*((ct(i-1)-c(i-1))/lambda-c(i-1))   ...
    -v3*c(i-1)^n3/(k3^n3+c(i-1)^n3) ) + jpm);
    ct(i)=ct(i-1)+dt*jpm;
    
end

atrans = find(t>=trans); t=t(atrans); c=c(atrans); ct=ct(atrans);
out.t=t; out.c=c; out.ct=ct;

[pks,locs] = findpeaks(ct);
ct_iei = diff(locs)*dt; out.ct_iei=ct_iei; out.pks = pks; out.locs = locs*dt+trans;
%ct_iei = diff(ct_cross)*dt; ct_iei = ct_iei(find(ct_iei>ct_ieimin)); out.ct_iei=ct_iei;

%% output & post analysis

mmc = [max(c) min(c)];

return

%% additional functions

function y = xinf(x,theta,k)
y=1./(1+exp(4*(theta-x)/k));
return

function f = finf(x,theta1,k1,theta2,k2)
f=1./(1+exp((theta1-x)/k1))./(1+exp((theta2-x)/k2));
return

function str = param_str(varargin)
gap0 = '+';
str =[];
for i=1:nargin
    if i==1, gap=''; else gap=gap0; end
    param_str = inputname(i);
    param_val = varargin{i};
    str = [ str gap param_str '=' num2str(param_val) ];
end

return

function [eupnea_i,sigh_i] = ClassifyEvents(pks,locs,ct)

    % Assume all bursts are inspiratory. Unless the ct has a negative
    % slope.
    
    Nvents = length(pks);
    
    eupnea_i = [];
    sigh_i   = [];
    
    for k = 1:Nvents
        
        event_i = locs(k);
        
        if ct(event_i-1) > ct(event_i+1) % then ct is going down
            sigh_i(end+1) = k; %#ok<*AGROW>
        else
            eupnea_i(end+1) = k;
        end
        
    end
    
    
    
    
return
















