%% Manuscript Figure #02 - Eupnea and Sigh Voltage Dependence
% script to make and export figures for figure 1 of manuscript
clear; close all; clc;

% figure stuff

%% first row, simulation, on eupnea system
includec = 0;
total = 3000;
fig = zeros(1,11); fig(2)=0;

thetaa_list = [-0.4:0.02:0]*10; % 0 is default 
[mean_a_iei,std_a_iei,mean_ct_iei,std_ct_iei,delta_s,delta_theta,delta_a]=deal(0*thetaa_list);

parfor i=1:length(thetaa_list)

thetaa=thetaa_list(i);    
    
    % IEI study
    [param, out] = tabakrinzelcalcium('thetaa',thetaa,...
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
    
    disp(['Finished run ' mat2str(i) ' out of ' mat2str(length(thetaa_list)) '.'])
                                    
end
                                 
%% Summary plot

%plot(thetaa_list,mean_a_iei,'o'); hold on;
%errorbar(thetaa_list,mean_a_iei,std_a_iei)
%plot(thetaa_list,mean_ct_iei,'o'); hold on;
%errorbar(thetaa_list,mean_ct_iei,std_ct_iei)


freq_e_iei=1./mean_e_iei;
freq_s_iei=1./mean_s_iei;

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

%% Export data

data_out = [thetaa_list' freq_e_iei' freq_s_iei'];
save('tikz/data/eup_sigh_freq_vs_thetaa_simulated.dat','data_out','-ascii')

function y = xinf(x,theta,k,lambda)
    y=lambda./(1+exp(4*(theta-x)/k));
end

function f = finf(x,theta1,k1,theta2,k2)
    f=1./(1+exp((theta1-x)/k1))./(1+exp((theta2-x)/k2));
end
