close all; clear; clc;

% load data
d = readtable('PKeupnea_sigh_freq.xlsx');

% establish concentrations
conc = 3:9;

% format data from each experiment and the means at each conc.
sighf = table2array(d(1:7,2:end));
msighf = mean(sighf,2,'omitnan');

eupnf = table2array(d(23:29,2:end));
meupnf = mean(eupnf,2,'omitnan');

% trying the other set of data, the "all bursts (burstlett) data"
hmf = table2array(d(12:18,2:end));
mhmf = mean(hmf,2,'omitnan');

% Further formatting
    % jitter
    jit = 0.00;
    A = repmat(conc,size(eupnf,2),1)';
    A = A + jit * randn([size(A)]);
    % shifts
    offset = 0.1;
    Aeup = A-offset;
    Asig = A+offset;
    conceup = conc-offset;
    concsig = conc+offset;
    
% Linear Regression
FOsigh = fit(conc',msighf,'poly1');
FOeupn = fit(conc',meupnf,'poly1');

% MATLAB plots
plot(concsig,msighf,'b--',conceup,meupnf,'r--',conc,mhmf,'m--'); hold on;
plot(Asig,sighf,'bx',Aeup,eupnf,'rx')
set(gca,'XDir','reverse')
legend('sigh','eup')

% Readying data for export to .dat
Bsig = Asig'; Beup = Aeup';
dataout_means = [concsig',conceup', msighf,meupnf];
dataout_eup = [Beup eupnf'];
dataout_sig = [Bsig sighf'];
% ready fit oject for export
FOsighL = FOsigh(conc);
FOeupnL = FOeupn(conc);
FOsighO = [conc' FOsighL];
FOeupnO = [conc' FOeupnL];

%dataout_eup(isnan(dataout_eup)) = 0;

% save data as .dat
save('tikz/data/PKdata_means.dat','dataout_means','-ascii')
save('tikz/data/PKdata_eup.dat','dataout_eup','-ascii')
save('tikz/data/PKdata_sig.dat','dataout_sig','-ascii')
save('tikz/data/PKdata_eupfit.dat','FOeupnO','-ascii')
save('tikz/data/PKdata_sigfit.dat','FOsighO','-ascii')