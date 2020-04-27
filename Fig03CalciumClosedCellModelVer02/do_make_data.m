%% Tabak-Rinzel-like model
%% do_make_data.m example

clc; clear; close all;
system('rm -r fig*')


%% Example calcium oscillation open cell model
[param, out] = tabakrinzelcalcium('writetraj',1,'filenametraj','figdata1',...
    'includec',1,'includetheta',0,'includes',0,...
    'scale_calcium',0.5,...
    'jin0',0.02,'jin1',0,...
    'trans',2000,'total',2150,'thintraj',100,'fig',[1 0 0 0 1 0 0 0 0 ],...
    'seed',-1);

% write ctnullcline .dat files
fid = fopen('fig_ctnull.dat','w');
fprintf(fid,'%6.4f  %6.4f  \n',[ out.ctnull_c ; out.ctnull_ct ]);
fclose(fid);

fid = fopen('fig_cnull.dat','w');
fprintf(fid,'%6.4f  %6.4f  \n',[ out.cnull_c ; out.cnull_ct ]);
fclose(fid);

v1=param.v1
v2=param.v2
v3=param.v3
k3=param.k3
snapnow

