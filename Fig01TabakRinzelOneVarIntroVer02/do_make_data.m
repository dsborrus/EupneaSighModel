%% Tabak-Rinzel-like model  
%% do_make_data.m example
 

% greg test change for github

clc; clear; close all;
system('rm -r fig*')
rng(201222)

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
sr = 10; %sampling rate

    
%% Example eupnea simulation with 1 slow var (type 2)
    [param, out] = tabakrinzelcalcium('writetraj',0,'filenametraj','figdata1',...
    'includec',0,'includetheta',0,...
    'scale_eupnea',scale_eupnea2,'N',N,......
    'trans',2000,'total',2025,'thintraj',10,'fig',[1 1 0 0 0 0 0 0 0 ],...
    'seed',0);

figdata1 = [out.t(1:sr:end)'-2000 out.a(1:sr:end)' out.s(1:sr:end)'];
save(['figs_and_data_' filename '/figdata1.dat'],'figdata1','-ascii')

snapnow

%% Example eupnea simulation with 1 slow var (type 1)
    [param, out] = tabakrinzelcalcium('writetraj',0,'filenametraj','figdata2',...
    'includec',0,'includetheta',0,...
    'scale_eupnea',scale_eupnea1,...
    'thetas',thetas1,'ks',ks1,'N',N,...
    'trans',2000,'total',2025,'thintraj',10,'fig',[1 1 0 0 0 0 0 0 0 ],...
    'seed',0,...
    'tausmax',tausmax1,'tausmin',tausmin1);

figdata2 = [out.t(1:sr:end)'-2000 out.a(1:sr:end)' out.s(1:sr:end)'];
save(['figs_and_data_' filename '/figdata2.dat'],'figdata2','-ascii')

snapnow


%% Histogram (type 2)
    [param, out] = tabakrinzelcalcium('writehist',0,'filenamehist','figdata3',...
    'includec',0,'includetheta',0,...
    'scale_eupnea',scale_eupnea2,'N',N,......
    'trans',1000,'total',5000,'fig',[0 0 0 0 0 0 1 0 0 ],...
    'seed',0);

out.mean_a_iei;
out.cv_a_iei;

% find inter burst intervals and cooresponding peak amplitude
[pks,locs] = findpeaks(out.a,'minpeakprominence',0.5);
iei = diff(out.t(locs));
pks2 = pks(2:end);

% make the histogram
h=histogram(iei);
BE = h.BinEdges(1:end-1)';
V  = h.Values';
figdata3hist = [BE V/sum(V)];
figdata3 = [iei' pks2'];

% export both
save(['figs_and_data_' filename '/figdata3.dat'],'figdata3','-ascii')
save(['figs_and_data_' filename '/figdata3hist.dat'],'figdata3hist','-ascii')

% make a fit of the histogram
gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
f1 = fit(BE,V/sum(V),'gauss1','start',[1 4 1]);
gaussfit = [(3.4:0.01:4.8)'+((BE(2)-BE(1))/2) f1(3.4:0.01:4.8)];%+((BE(2)-BE(1))/2) shifts the fit a bit to the right, to make it centered with the bars
save(['figs_and_data_' filename '/figdata3gaussfit.dat'],'gaussfit','-ascii')
% make a fit of the iei data against burst size
f1 = fit(iei',pks2','poly1');
polyfit = [(3.4:0.1:4.8)' f1(3.4:0.1:4.8)];
save(['figs_and_data_' filename '/figdata3polyfit.dat'],'polyfit','-ascii')

snapnow

%% Histogram (type 1)
    [param, out] = tabakrinzelcalcium('writehist',0,'filenamehist','figdata4',...
    'includec',0,'includetheta',0,...
    'scale_eupnea',scale_eupnea1,...
    'thetas',thetas1,'ks',ks1,'N',N,...
    'trans',1000,'total',5000,'fig',[0 0 0 0 0 0 1 0 0 ],...
    'seed',0,...
    'tausmax',tausmax1,'tausmin',tausmin1);
out.mean_a_iei;
out.cv_a_iei;
% find inter burst intervals and cooresponding peak amplitude
[pks,locs] = findpeaks(out.a,'minpeakprominence',0.5);
iei = diff(out.t(locs));
pks2 = pks(2:end);

% make the histogram
h=histogram(iei);
BE = h.BinEdges(1:end-1)';
V  = h.Values';
figdata4hist = [BE V/sum(V)];
figdata4 = [iei' pks2'];

% export both
save(['figs_and_data_' filename '/figdata4.dat'],'figdata4','-ascii')
save(['figs_and_data_' filename '/figdata4hist.dat'],'figdata4hist','-ascii')

% make a fit of the histogram
expoEqn = 'A*exp(r*x)';
f1 = fit(BE,V/sum(V),expoEqn,'start',[1 -1]);
gaussfit = [(-1:0.1:30)'+((BE(2)-BE(1))/2) f1(-1:0.1:30)]; %+((BE(2)-BE(1))/2) shifts the fit a bit to the right, to make it centered with the bars
save(['figs_and_data_' filename '/figdata4gaussfit.dat'],'gaussfit','-ascii')
% make a fit of the iei data against burst size
f1 = fit(iei',pks2','poly1');
polyfit = [(0:5:30)' f1(0:5:30)];
save(['figs_and_data_' filename '/figdata4polyfit.dat'],'polyfit','-ascii')

snapnow


%% Two Pulse Protocol (type 2)
    [param, out] = tabakrinzelcalcium('writetraj',0,'filenametraj','figdata5',...
    'includec',0,'includetheta',0,...
    'scale_eupnea',scale_eupnea2,'N',N,......
    'trans',0,'total',7,'thintraj',10,'fig',[1 1 0 0 0 0 0 0 0 ],...
    'seed',0,...
    'init_a',0.3,'init_s',0.8,'jump_a',0.2,'jump_t',[5.5]);

figdata5 = [out.t(1:sr:end)' out.a(1:sr:end)' out.s(1:sr:end)'];
save(['figs_and_data_' filename '/figdata5.dat'],'figdata5','-ascii')

snapnow

%% Two Pulse Protocol (type 1)
  [param, out] = tabakrinzelcalcium('writetraj',0,'filenametraj','figdata6',...
    'includec',0,'includetheta',0,...
    'scale_eupnea',scale_eupnea1,...
    'thetas',thetas1,'ks',ks1,'N',N,...
    'trans',0,'total',7,'thintraj',10,'fig',[1 1 0 0 0 0 0 0 0 ],...
    'seed',3,...
    'init_a',0.3,'init_s',0.8,'jump_a',0.2,'jump_t',[5],...
    'tausmax',tausmax1,'tausmin',tausmin1);

figdata6 = [out.t(1:sr:end)' out.a(1:sr:end)' out.s(1:sr:end)'];
save(['figs_and_data_' filename '/figdata6.dat'],'figdata6','-ascii')

snapnow

%% clean up, save workspace

save(['figs_and_data_' filename '/workspace.mat'])