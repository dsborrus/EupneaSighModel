%% Tabak-Rinzel-like model
%% do_make_data.m 
% Interested in the effect of calcium dynanmics on event coupling 
% this is the script to actually make the fig we might use in the paper

% version 1 - everythin assembled. Stopped when sigh/eupnea classification
% changed. In this version it is amplitude (200408)

% version 2 - opened 200408. Modifying sigh/eupnea classifcation. Changed a
% to go from 0 to amax, instead of 0 to 1.

clc; clear; close all; mytic=tic;

name = ['PubFig_' char(datetime('now','Format','yyMMdd_hhmmss'))];
system(['mkdir ' name]);
diary([name '/diaryfile'])
disp(['Running publiciation figure: ' name])

diary([name '/diaryfile'])

% % parameters % %

tmax = 30000;
includec = 1;
scale_calcium = 0.5;
dosave = 0;
N = 80;

amax=25;

%% control conditions
jin0 = scale_calcium*0.008; % 0.02 0.018 0
jin1=scale_calcium/amax * 0.2; % 0.02 0.04 0.2

[paramsC, outC] = tabakrinzelcalcium('ESCOUPLINGvsCa',1,...
                                   'ESCouplingFig',0,...
                                   'total',tmax,...
                                   'includec',includec,...
                                   'N',N,...
                                   'jin0',jin0,...
                                   'jin1',jin1,...
                                   'fig',[0 0 0 0 0 0 0 0 0],...
                                   'dosave',dosave...
                                   );

%% larger jin1 conditions (1.5x larger)
jin1D1=jin1*(1.5); % 0.02 0.04 0.2

[paramsD1, outD1] = tabakrinzelcalcium('ESCOUPLINGvsCa',1,...
                                       'ESCouplingFig',0,...
                                       'total',tmax,...
                                       'includec',includec,...
                                       'N',N,...
                                       'jin0',jin0,...
                                       'jin1',jin1D1,...
                                       'fig',[0 0 0 0 0 0 0 0 0],...
                                       'dosave',dosave...
                                       );
                                   
%% smaller jin1 conditions (0.5x the size)
jin1D2=jin1*(0.5); % 0.02 0.04 0.2

[paramsD2, outD2] = tabakrinzelcalcium('ESCOUPLINGvsCa',1,...
                                       'ESCouplingFig',0,...
                                       'total',tmax,...
                                       'includec',includec,...
                                       'N',N,...
                                       'jin0',jin0,...
                                       'jin1',jin1D2,...
                                       'fig',[0 0 0 0 0 0 0 0 0],...
                                       'dosave',dosave...
                                       );                                   
                               
%% figure creation                               
                               
figure('position',[1000 100 1500 1000])

twin = 200;

% control plot
subplot(3,3,2); hold on;
plot(outC.t,outC.a,'k');
plot(outC.t(outC.Coupling.locs(outC.Coupling.sigh_i)),outC.Coupling.pks(outC.Coupling.sigh_i),'ro')
plot(outC.t(outC.Coupling.locs(outC.Coupling.sigh_i-1)),outC.Coupling.pks(outC.Coupling.sigh_i-1),'bo')
xlim([2000 2000+twin])
subplot(3,3,5);                        
histogram(outC.Coupling.EuSiIn)

% larger jin1 plot
subplot(3,3,3); hold on;
plot(outD1.t,outD1.a,'k');
plot(outD1.t(outD1.Coupling.locs(outD1.Coupling.sigh_i)),outD1.Coupling.pks(outD1.Coupling.sigh_i),'ro')
plot(outD1.t(outD1.Coupling.locs(outD1.Coupling.sigh_i-1)),outD1.Coupling.pks(outD1.Coupling.sigh_i-1),'bo')
xlim([2000 2000+twin])
subplot(3,3,6);                        
histogram(outD1.Coupling.EuSiIn)
                               
% smaller jin1 plot
subplot(3,3,1); hold on;
plot(outD2.t,outD2.a,'k');
plot(outD2.t(outD2.Coupling.locs(outD2.Coupling.sigh_i)),outD2.Coupling.pks(outD2.Coupling.sigh_i),'ro')
plot(outD2.t(outD2.Coupling.locs(outD2.Coupling.sigh_i-1)),outD2.Coupling.pks(outD2.Coupling.sigh_i-1),'bo')
xlim([2000 2000+twin])
subplot(3,3,4);                        
histogram(outD2.Coupling.EuSiIn)
%% Clean up 
%save([sweepname '/matfile.mat'])
t_taken = toc(mytic);
%if ispc; system(['move *.png ' name '/']); end   
disp(['Time taken was ~' mat2str(round(t_taken,0)) ' seconds. Or ~' mat2str(round(t_taken/60,0)) ' minute.' ])
disp('Exiting script. Have a nice day.')
diary off

%% Extra functions