close all; clc; clear;

load('exampleeupnea.mat');
t_pre2 = readtable('exampleeupnea.txt');
t_pre = table2array(t_pre2(:,1));

ds = 1/1000; % downsample

a = round(length(C3B3)*ds);

y = C3B3(1:a:end);
t = t_pre(1:a:end);

out = [t y];

%plot(t,y)
save('sampletrace.dat','out','-ascii')