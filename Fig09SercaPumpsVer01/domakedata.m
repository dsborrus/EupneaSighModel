% Script to run model with different blocks of serca
% We'll need three cases...
% no block, partial block, total block

% This work mirrors the thapsigargin experiments.
% This script will get the data and export it for latex

% begin script
clear; close all;

% consistent parameters
tmax=2200;trans=2000;
figs = [1 0 0 0 0 0 0 0 0];
dt = 0.001;

% changing parameter
% serca pump (v3) = [100% 60% 0%]
v3 = [60 42 12]

% initialize some things
n = length(v3);
t = zeros((tmax-trans)/dt+1,1);
a = zeros((tmax-trans)/dt+1,3); ct=a; c=a;

% loop through
for k = 1:n
    [p,o] = tabakrinzelcalcium('total',tmax,...
                               'trans',trans,...
                               'fig',figs,...
                               'dt',dt,...
                               'v3',v3(k));
    
    a(:,k)  = o.a;
    ct(:,k) = o.ct;
    c(:,k)  = o.c;
end
t = o.t-trans;

figure
plot(t,a(:,1)); hold on;
plot(t,a(:,2))
plot(t,a(:,3))


%system('rm -r fig*')