%% Manuscript Figure #02 - Eupnea and Sigh Voltage Dependence
% script to make and export figures for figure 1 of manuscript
clear; close all; clc;

% figure stuff

%% first row, simulation, on eupnea system
total = 3000; % total time of simulation
jin0 = 0.005; % default 0.0045 0.0036
jin1 = 0.00512; % default 0.0004 0.00512


fig = zeros(1,11); fig(2)=0;

thetaa_list = [-0.3:0.04:0]*10; % 0 is default 
[mean_a_iei,std_a_iei,mean_ct_iei,std_ct_iei,delta_s,delta_theta,delta_a]=deal(0*thetaa_list);

% special simulations to grab
% how about data from thetaa = -3 and thetaa=0. Or the min and max
downsampleN = 100;
maxwindow = 240/0.001;

for i=1:length(thetaa_list)

    thetaa=thetaa_list(i);    
    
    % IEI study
    [param, out] = tabakrinzelcalcium('thetaa',thetaa,...
                                      'jin0',jin0,'jin1',jin1,...
                                      'tauthetamax',8,...
                                      'trans',1000,...
                                      'total',total,...
                                      'fig',[0 0 0 0 0 0 0 0 0 ],...
                                      'ESCOUPLINGvsCa',1 ...
                                      );

    
    % individual figures (non default)
    if 0
        figure
        plot(out.t,out.a,'b')
        hold on;
        plot(out.t(out.Coupling.locs(out.Coupling.eupnea_i)),5,'ro');
        plot(out.t(out.Coupling.locs(out.Coupling.sigh_i)),7,'ko');
    end
    
    mean_e_iei(i) = mean(diff(out.t(out.Coupling.locs(out.Coupling.eupnea_i)))); 
    std_e_iei(i)  = std(diff(out.t(out.Coupling.locs(out.Coupling.eupnea_i))));
    mean_s_iei(i) = mean(diff(out.t(out.Coupling.locs(out.Coupling.sigh_i)))); 
    std_s_iei(i)  = std(diff(out.t(out.Coupling.locs(out.Coupling.sigh_i))));
    
    if i == 1
        excited_simulation = [out.t(1:downsampleN:maxwindow)'-out.t(1) out.a(1:downsampleN:maxwindow)'];
        save('tikz/data/excited_simulation.dat','excited_simulation','-ascii')
        disp(['excited simulations thetaa is ' mat2str(param.thetaa)])
    elseif i == length(thetaa_list)
        depressed_simulation = [out.t(1:downsampleN:maxwindow)'-out.t(1) out.a(1:downsampleN:maxwindow)'];
        save('tikz/data/depressed_simulation.dat','depressed_simulation','-ascii')
        disp(['depressed simulations thetaa is ' mat2str(param.thetaa)])
    end
    
    disp(['Finished run ' mat2str(i) ' out of ' mat2str(length(thetaa_list)) '.'])
    
    
                                    
end

% IEI study
[param, out] = tabakrinzelcalcium('thetaa',thetaa,...
                                  'jin0',jin0,'jin1',jin1,...
                                  'tauthetamax',8,...
                                  'trans',1000,...
                                  'total',total,...
                                  'fig',[0 0 0 0 1 0 0 0 0 ],...
                                  'ESCOUPLINGvsCa',1 ...
                                  );
                              
disp('Running once more for the nulc lines')


% write ctnullcline .dat files
fid = fopen('tikz/data/fig_cnull.dat','w');
fprintf(fid,'%6.4f  %6.4f  \n',[ out.ctnull_c ; out.ctnull_ct ]);
fclose(fid);

fid = fopen('tikz/data/fig_ctnull_low.dat','w');
fprintf(fid,'%6.4f  %6.4f  \n',[ out.cnull_c ; out.cnull_ct ]);
fclose(fid);
fid = fopen('tikz/data/fig_ctnull_high.dat','w');
fprintf(fid,'%6.4f  %6.4f  \n',[ out.cnull_c+0.05 ; out.cnull_ct ]);
fclose(fid);

ctnull_low = ((param.k4*(param.jin0+param.jin1*0))/(param.v4-param.jin0+param.jin1*0))^(1/4);
ctnull_high = ((param.k4*(param.jin0+param.jin1*10))/(param.v4-param.jin0+param.jin1*10))^(1/4);


                                 
%% Summary plot

%plot(thetaa_list,mean_a_iei,'o'); hold on;
%errorbar(thetaa_list,mean_a_iei,std_a_iei)
%plot(thetaa_list,mean_ct_iei,'o'); hold on;
%errorbar(thetaa_list,mean_ct_iei,std_ct_iei)


freq_e_iei=1./mean_e_iei;
freq_s_iei=60./mean_s_iei;

figure
subplot(2,1,1)
plot(thetaa_list,freq_e_iei,'o'); hold on;
xlabel('\theta_a')
ylabel('eupnea freq (Hz)')

subplot(2,1,2)
plot(thetaa_list,freq_s_iei,'o'); hold on;
xlabel('\theta_a')
ylabel('sigh freq (Hz)')


figure
plot(thetaa_list,freq_e_iei,'o'); hold on;
plot(thetaa_list,freq_s_iei,'o');
xlabel('\theta_a')
ylabel('freq (Hz)')
legend('Eupnea','Sigh')

disp(['Max eupnea frequency = ' mat2str(max(freq_e_iei))])
disp(['Min eupnea frequency = ' mat2str(min(freq_e_iei))])
disp(['Max sigh frequency = ' mat2str(max(freq_s_iei))])
disp(['Min sigh frequency = ' mat2str(min(freq_s_iei))])

figure
curtime = datetime('now','format','yyyyMMdd-HHmmss');
plot(thetaa_list,freq_e_iei/max(freq_e_iei),'bo','linewidth',2); hold on;
plot(thetaa_list,freq_s_iei/max(freq_s_iei),'ro','linewidth',2);
xlabel('\theta_a')
ylabel('freq (relative to maximum)')
title(['j0 = ' mat2str(jin0) '.   j1 = ' mat2str(jin1) '.   File: ' char(curtime)])
legend('Eupnea','Sigh')
saveas(gcf,['ParamSweepData/Simulation-' char(curtime) '.png'])

%% Export data

data_out = [thetaa_list' freq_e_iei' freq_s_iei'];
%save(['tikz/data/eup_sigh_freq_vs_thetaa_simulated_v' char(curtime) '.dat'],'data_out','-ascii')
save('tikz/data/eup_sigh_freq_vs_thetaa_simulated.dat','data_out','-ascii')

function y = xinf(x,theta,k,lambda)
    y=lambda./(1+exp(4*(theta-x)/k));
end

function f = finf(x,theta1,k1,theta2,k2)
    f=1./(1+exp((theta1-x)/k1))./(1+exp((theta2-x)/k2));
end
