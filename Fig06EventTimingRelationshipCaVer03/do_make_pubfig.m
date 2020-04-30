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

tmax = 5000;
includec = 1;
scale_calcium = 0.5;
dosave = 0;
N = 80;

amax=25;

% control conditions
jin0 = scale_calcium*0.019; % 0.02 0.018 0
jin1=scale_calcium/amax * 0.1; % 0.02 0.04 0.2

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
                               

                               
figure('position',[100 100 800 600])                               
subplot(3,3,5);                               
                               
%% Clean up 
%save([sweepname '/matfile.mat'])
t_taken = toc(mytic);
if ispc; system(['move *.png ' name '/']); end   
disp(['Time taken was ~' mat2str(round(t_taken,0)) ' seconds. Or ~' mat2str(round(t_taken/60,0)) ' minute.' ])
disp('Exiting script. Have a nice day.')
diary off

%% Extra functions