%% Tabak-Rinzel-like model
%% do_make_data.m 
% Exploring efffect of excitability 
% on the two-variable version of he Tabak-Rinzel-like model. 
% Here \theta is a dynamic variable.
% Changes in network excitability modeled by changing \thetaa. 
% Looking at IEIs of eupnea and sigh 

clc; clear; close all;
system('rm -r fig*')

thetaa_list = [ -0.4:0.02:0.1 ]; % 0 is default 
[mean_a_iei,std_a_iei,mean_ct_iei,std_ct_iei,delta_s,delta_theta,delta_a]=deal(0*thetaa_list);
 for i=1:length(thetaa_list)
    
    thetaa=thetaa_list(i);
    
    % Example phase plane
    [param, out] = tabakrinzelcalcium('thetaa',thetaa,'trans',1000,'total',1100,'fig',[1 0 0 1 0 0 0 1 1]);
    snapnow
    
    % IEI study
    [param, out] = tabakrinzelcalcium('thetaa',thetaa,'trans',1000,'total',2000,'fig',[0 0 0 0 0 0 1 0 0 ]);
    snapnow
    
    mean_a_iei(i)=out.mean_a_iei; std_a_iei(i)=out.std_a_iei;
    mean_ct_iei(i)=out.mean_ct_iei; std_ct_iei(i)=out.std_ct_iei;
    
    %delta_s(i)=out.delta_s;
    %delta_theta(i)=out.delta_theta;
    %delta_a(i)=out.delta_a;
end



%% Summary plot

%plot(thetaa_list,mean_a_iei,'o'); hold on;
%errorbar(thetaa_list,mean_a_iei,std_a_iei)
%plot(thetaa_list,mean_ct_iei,'o'); hold on;
%errorbar(thetaa_list,mean_ct_iei,std_ct_iei)


freq_a_iei=1./mean_a_iei;
freq_ct_iei=1./mean_ct_iei;

figure(1)
subplot(2,1,1)
plot(thetaa_list,freq_a_iei,'o'); hold on;
xlabel('\theta_a')
ylabel('eupnea freq (Hz)')

subplot(2,1,2)
plot(thetaa_list,freq_ct_iei,'o'); hold on;
xlabel('\theta_a')
ylabel('sigh freq (Hz)')
snapnow

fid = fopen('freq_a_iei.dat','w');
fprintf(fid,'%6.4f  %6.4f  \n',[ thetaa_list; freq_a_iei ]);
fclose(fid);

fid = fopen('freq_ct_iei.dat','w');
fprintf(fid,'%6.4f  %6.4f  \n',[ thetaa_list; freq_ct_iei ]);
fclose(fid);


return


%% Clean up 
save('html/param_study.mat')
system(['mv fig* html' ]); % collect runs (each is a directory) 

%% Version of TabakRinzelCalcium used 
type tabakrinzelcalcium

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

