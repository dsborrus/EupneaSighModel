%% Tabak-Rinzel-like model
%% do_make_data.m 
% Interested in the effect of calcium dynanmics on event coupling 
% start a parpool (ideally before script runs)
% this script runs the psweep for a range of jin0 and jin1.
% then plots a grid of figures for different values, including trajectories
% and distribution of coupling events.

% version 1 - everythin assembled. Stopped when sigh/eupnea classification
% changed. In this version it is amplitude (200408)

% version 2 - opened 200408. Modifying sigh/eupnea classifcation. Changed a
% to go from 0 to amax, instead of 0 to 1.

clc; clear; close all; tic; t_start=tic; 
% record
sweepname = ['Sweep_' char(datetime('now','Format','yyMMdd_hhmmss'))];
system(['mkdir ' sweepname]);
diary([sweepname '/diaryfile'])

% static nondefault parameters
tmax = 15000;
dt = 0.001;
includec = 1;
scale_calcium = 0.5;
compactout = 0;
N = 80;
fig=[0 0 0 0 0 0 0 0 0];

% Fraction of "activity" to be synaptic activity and fraction to be calcium
amax = 25; % note for user. amax here must be same as in TabakRinzelCalc model
% it can't be set on it's own just yet.

disp('param sweep range is:')
% base intake of Ca (0.18 at first) (Greg's notes 0.02 0.018 0)
jin0_min = 0.016 %#ok<*NOPTS>
jin0_max = 0.022
jin0_stp = 0.001
jin0_array = scale_calcium*(jin0_min:jin0_stp:jin0_max);
% Ca intake from eupnea (0.04 at first) (Greg's notes 0.02 0.04 0.2)
jin1_min = 0.03
jin1_max = 0.11
jin1_stp = 0.01
jin1_array = scale_calcium/amax*(jin1_min:jin1_stp:jin1_max);


% initialize
Nj0 = length(jin0_array);
Nj1 = length(jin1_array);
Ncases = Nj0*Nj1;
outs   = cell(Ncases,1);
params = cell(Ncases,1);
[swpX,swpY] = meshgrid(jin0_array,jin1_array);


% begin psweep
disp(['Beginning sweep. name=' sweepname])

parfor k = 1:Ncases
    
    jin0 = swpX(k);
    jin1 = swpY(k);

    [params{k}, outs{k}] = tabakrinzelcalcium('ESCOUPLINGvsCa',1,...
                                              'ESCouplingFig',0,...
                                              'includec',includec,...
                                              'total',tmax,...
                                              'dt',dt,...
                                              'N',N,...
                                              'jin0',jin0,...
                                              'jin1',jin1,...
                                              'compactout',compactout,...
                                              'fig',fig);
    
    %waitbar stuff                                                      
    %send(D, k);
    %waitbar(k/Ncases,h)
    disp(['Simulation #' mat2str(k) ' done of ' mat2str(Ncases) ' total.'])
                                                          
end
    

%% Summary plot
close all
% actual activity grid plots
if 1 && compactout==0
    disp('Starting figure creation. Making trajectory grid figure')
    F1 = figure('units','normalized','outerposition',[0 0 1 1]);
    trjwnd = 200;
    %[fX,fY] = meshgrid(1:Nj0,1:Nj1); 
    for fi = 1:Ncases
        subplot(Nj0,Nj1,fi); hold on;
        
        simots = outs{fi};
        simprms = params{fi};
        
        t=simots.t;
        a=simots.a;
        locs = simots.Coupling.locs;
        pks = simots.Coupling.pks;
        eupnea_i = simots.Coupling.eupnea_i;
        sigh_i = simots.Coupling.sigh_i;
        
        eupnea_i = eupnea_i(find(locs(eupnea_i)<(trjwnd*1000)));
        sigh_i = sigh_i(find(locs(sigh_i)<(trjwnd*1000)));
                
        plot(t(1:trjwnd*1000),a(1:trjwnd*1000),'b')
        plot(t(locs(sigh_i)),pks(sigh_i),'ro')
        plot(t(locs(sigh_i-1)),pks(sigh_i-1),'go')
        title(['jin0 = ' mat2str(simprms.jin0/simprms.scale_calcium) ', jin1 = ' mat2str(simprms.jin1*amax/simprms.scale_calcium)])
        
    end
    
    savefig(F1,[sweepname '/Trajectories'])
    saveas(F1,[sweepname '/Trajectories.png'])
    
end

if 1
    disp('Starting figure creation. Making histogram grid figure')
    F2 = figure('units','normalized','outerposition',[0 0 1 1]);
    %[fX,fY] = meshgrid(1:Nj0,1:Nj1); 
    for fi = 1:Ncases
        subplot(Nj0,Nj1,fi); hold on;
        
        simots = outs{fi};
        simprms = params{fi};
        
        histogram(simots.Coupling.EuSiIn,'binedges',[0:0.05:1.5]);
        title(['jin0 = ' mat2str(simprms.jin0/simprms.scale_calcium) ', jin1 = ' mat2str(simprms.jin1*amax/simprms.scale_calcium)])
        if fi==1; xlabel('Inspiratry-Sigh interval'); end
        
        % annotate graph w/ number of eupneas
        if 1
            mNepn = mean(diff(simots.Coupling.sigh_i));
            sNepn = std(diff(simots.Coupling.sigh_i));
            dim = [.1 .1 .3 .3];
            str = ['N eup = ' mat2str(round(mNepn,2)) ' pm ' mat2str(round(sNepn,2))];
            text(0.6,0.9,str,'units','normalized')
            %annotation('textbox','String',str); 
        end
        
    end
    savefig(F2,[sweepname '/HistogramGrids'])
    saveas(F2,[sweepname '/HistogramGrids.png'])
end
toc
%% Clean up 
%save([sweepname '/matfile.mat'])
t_taken = toc(t_start);
disp(['Parameter sweep complete. Output saved to ' sweepname '.'])
disp(['Time taken was ~' mat2str(round(t_taken,0)) ' seconds. Or ~' mat2str(round(t_taken/60,0)) ' minute.' ])

disp('Exiting script. Have a nice day.')
diary off
%% Extra functions