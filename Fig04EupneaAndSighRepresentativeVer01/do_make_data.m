%% Tabak-Rinzel-like model
%% do_make_data.m example

clc; clear; close all;
system('rm -r fig*')

jin0=0;
jin1=0.1;
thetac=0.33;

%% Example calcium oscillation open cell model
[param, out] = tabakrinzelcalcium('writetraj',1,'filenametraj','figdata1',...
    'includec',1,'includetheta',1,'includes',1,...
    'trans',3000,'total',3300,'thintraj',100,'fig',[1 0 0 0 1 0 0 0 0 ],...
    'seed',0,'jin0',jin0,'jin1',jin1,'thetac',thetac);

[param, out] = tabakrinzelcalcium('writetraj',1,'filenametraj','figdata2',...
    'includec',1,'includetheta',1,'includes',1,...
    'trans',3120,'total',3160,'thintraj',10,'fig',[1 0 0 0 1 0 0 0 0 ],...
    'seed',0,'jin0',jin0,'jin1',jin1,'thetac',thetac);

return


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

