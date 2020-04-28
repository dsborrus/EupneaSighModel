%% Tabak-Rinzel-like model

clc; clear; close all;

tmax = 2200;
scale_calcium = 0.5;
dosave = 0;
NoiseTuning_rhythmfreqfig = 1;
N = 50000;

amax=12.5;

jin0 = scale_calcium*0.018; % 0.02 0.018 0
jin1=scale_calcium/amax * 0.04; % 0.02 0.04 0.2

name = ['OneOff_' char(datetime('now','Format','yyMMdd_hhmmss'))];
system(['mkdir ' name]);
diary([name '/diaryfile'])
disp(['Running one off: ' name])


[params, out] = tabakrinzelcalcium('ESCOUPLINGvsCa',0,...
                                   'ESCouplingFig',0,...
                                   'total',tmax,...
                                   'N',N,...
                                   'jin0',jin0,...
                                   'jin1',jin1,...
                                   'fig',[1 0 0 0 0 0 0 0 0],...
                                   'dosave',dosave,...
                                   'NoiseTuning_rhythmfreqfig',NoiseTuning_rhythmfreqfig...
                                   );
            
    
                               

if dosave                               
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = get(FigHandle, 'Name');
      savefig(FigList(iFig), [name '/Fig' mat2str(iFig)] );
    end 
end

diary off