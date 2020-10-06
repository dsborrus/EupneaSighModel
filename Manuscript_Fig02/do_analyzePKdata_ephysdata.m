close all; clear; clc;

d = readtable('PKeupnea_sigh_freq.xlsx');

conc = 3:9;

sighf = table2array(d(1:7,2:end))/60;
msighf = mean(sighf,2,'omitnan');

eupnf = table2array(d(23:29,2:end));
meupnf = mean(eupnf,2,'omitnan');

% jitter
jit = 0.05;
A = repmat(conc,size(eupnf,2),1)';
A = A + jit * randn([size(A)]);

plot(conc,msighf,'b--',conc,meupnf,'r--'); hold on;
plot(A,sighf,'bx',A,eupnf,'rx')
set(gca,'XDir','reverse')

B = A';
dataout_means = [conc', msighf,meupnf];
dataout_eup = [B eupnf'];
dataout_sig = [B sighf'];

%dataout_eup(isnan(dataout_eup)) = 0;

save('tikz/data/PKdata_means.dat','dataout_means','-ascii')
save('tikz/data/PKdata_eup.dat','dataout_eup','-ascii')
save('tikz/data/PKdata_sig.dat','dataout_sig','-ascii')
