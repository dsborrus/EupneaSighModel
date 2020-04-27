%% Tabak-Rinzel-like model  
%% do_make_data.m example
 

% greg test change for github

clc; clear; close all;
system('rm -r fig*')

tausmax1=0.1;
tausmin1=0.1;

scale_eupnea2=0.6;
scale_eupnea1=0.3;
N=100;
thetas1=0.3;
ks1=-1.0;
    
%% Example eupnea simulation with 1 slow var (type 2)
    [param, out] = tabakrinzelcalcium('writetraj',1,'filenametraj','figdata1',...
    'includec',0,'includetheta',0,...
    'scale_eupnea',scale_eupnea2,'N',N,......
    'trans',2000,'total',2025,'thintraj',10,'fig',[1 1 0 0 0 0 0 0 0 ],...
    'seed',0);
snapnow

%% Example eupnea simulation with 1 slow var (type 1)
    [param, out] = tabakrinzelcalcium('writetraj',1,'filenametraj','figdata2',...
    'includec',0,'includetheta',0,...
    'scale_eupnea',scale_eupnea1,...
    'thetas',thetas1,'ks',ks1,'N',N,...
    'trans',2000,'total',2025,'thintraj',10,'fig',[1 1 0 0 0 0 0 0 0 ],...
    'seed',0,...
    'tausmax',tausmax1,'tausmin',tausmin1);
snapnow


%% Histogram (type 2)
    [param, out] = tabakrinzelcalcium('writehist',1,'filenamehist','figdata3',...
    'includec',0,'includetheta',0,...
    'scale_eupnea',scale_eupnea2,'N',N,......
    'trans',1000,'total',5000,'fig',[0 0 0 0 0 0 1 0 0 ],...
    'seed',0);

out.mean_a_iei
out.cv_a_iei
snapnow

%% Histogram (type 1)
    [param, out] = tabakrinzelcalcium('writehist',1,'filenamehist','figdata4',...
    'includec',0,'includetheta',0,...
    'scale_eupnea',scale_eupnea1,...
    'thetas',thetas1,'ks',ks1,'N',N,...
    'trans',1000,'total',5000,'fig',[0 0 0 0 0 0 1 0 0 ],...
    'seed',0,...
    'tausmax',tausmax1,'tausmin',tausmin1);
out.mean_a_iei
out.cv_a_iei
snapnow


%% Two Pulse Protocol (type 2)
    [param, out] = tabakrinzelcalcium('writetraj',1,'filenametraj','figdata5',...
    'includec',0,'includetheta',0,...
    'scale_eupnea',scale_eupnea2,'N',N,......
    'trans',0,'total',7,'thintraj',10,'fig',[1 1 0 0 0 0 0 0 0 ],...
    'seed',0,...
    'init_a',0.3,'init_s',0.8,'jump_a',0.2,'jump_t',[5.5]);
snapnow

%% Two Pulse Protocol (type 1)
  [param, out] = tabakrinzelcalcium('writetraj',1,'filenametraj','figdata6',...
    'includec',0,'includetheta',0,...
    'scale_eupnea',scale_eupnea1,...
    'thetas',thetas1,'ks',ks1,'N',N,...
    'trans',0,'total',7,'thintraj',10,'fig',[1 1 0 0 0 0 0 0 0 ],...
    'seed',3,...
    'init_a',0.3,'init_s',0.8,'jump_a',0.2,'jump_t',[5],...
    'tausmax',tausmax1,'tausmin',tausmin1);
snapnow