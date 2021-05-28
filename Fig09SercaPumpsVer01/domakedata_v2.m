% Script to run model with different blocks of serca
% We'll need three cases...
% no block, partial block, total block

% This work mirrors the thapsigargin experiments.
% This script will get the data and export it for latex

% begin script
clear; close all;

% consistent parameters
tmax=3200;trans=2000;
figs = [0 0 0 0 0 0 0 0 0];
dt = 0.001;

% changing parameter
% serca pump (v3) = [100% 60% 0%]
v3array = [60 30 12]
t2c = [2400 2800];

[p,o] = tabakrinzelcalcium('total',tmax,...
                           'trans',trans,...
                           'fig',figs,...
                           'dt',dt,...
                           'v3array',v3array,...
                           't2c',t2c);

t = o.t-trans;

outdata = [t' o.a' o.c' o.ct'];

reduction=50;

outdatar = outdata(1:reduction:end,:);

save('outdata.dat','outdata','-ascii')
save('outdatar.dat','outdatar','-ascii')

if ismac
    system('cp *.dat ~/Library/Mobile\ Documents/com\~apple\~CloudDocs/My_Data/SighModelingStudy/Thapsigargin/ManuscriptFigure/Panel_model/');
end