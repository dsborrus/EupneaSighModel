%% Tabak-Rinzel-like model  

%% Tabak-Rinzel-like model  
%% do_make_data.m example
 

% greg test change for github

clc; clear; close all;
system('rm -r fig*')
rng(201222);

tausmax1=0.1;
tausmin1=0.1;

scale_eupnea2=0.6;
scale_eupnea1=0.3;
N=100;
thetas1=0.3;
ks1=-1.0;

% preparing directory for data
filename = char(datetime('now','format','yyMMdd_HHmmss'));
system(['mkdir figs_and_data_' filename]);

% fig params
sr = 100; %sampling rate

    
%% Example eupnea simulation with 2 slow var (type 2)
[param, out] = tabakrinzelcalcium('trans',2000,'total',2025,...
                                   'N',N,...
                                   'includec',0,'includetheta',1,...
                                   'fig',[1 1 0 0 0 0 0 0 0 ],...
                                   'seed',0);

figdata1b = [out.t(1:sr:end)'-2000 out.a(1:sr:end)' out.s(1:sr:end)' out.theta(1:sr:end)' zeros(length(out.a(1:sr:end)),1)];
save(['figs_and_data_' filename '/figdata1b.dat'],'figdata1b','-ascii')

%% Longer trace to measure inter event interval characteristics

[param, out] = tabakrinzelcalcium('trans',2000,'total',3000,...
                                   'N',N,...
                                   'includec',0,'includetheta',1,...
                                   'seed',0);

% find inter burst intervals and cooresponding peak amplitude
[pks,locs] = findpeaks(out.a,'minpeakprominence',0.5);
iei = diff(out.t(locs));
pks2 = pks(2:end);

% make the histogram
h=histogram(iei);
BE = h.BinEdges(1:end-1)';
V  = h.Values';
figdata2hist = [BE V/sum(V)];
figdata2 = [iei' pks2'];

% export both
save(['figs_and_data_' filename '/figdata2.dat'],'figdata2','-ascii')
save(['figs_and_data_' filename '/figdata2hist.dat'],'figdata2hist','-ascii')

% make a fit of the histogram
gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
f1 = fit(BE,V/sum(V),'gauss1','start',[1 4 1]);
gaussfit = [(4:0.01:5)'+((BE(2)-BE(1))/2) f1(4:0.01:5)]; %+((BE(2)-BE(1))/2) shifts the fit a bit to the right, to make it centered with the bars
save(['figs_and_data_' filename '/figdata2gaussfit.dat'],'gaussfit','-ascii')
% make a fit of the iei data against burst size
f1 = fit(iei',pks2','poly1');
polyfit = [(4:.1:5)' f1(4:.1:5)];
save(['figs_and_data_' filename '/figdata2polyfit.dat'],'polyfit','-ascii')                               

save(['figs_and_data_' filename '/workspace.mat'])