% this script runs the 2 parameter bifurcation diagram, but the poor man's
% bif diagram.

% We need to loop through all the different trajectories with different j0
% and j1. Then run a method to detect the appearence of eupnea and/or sigh.

% Version 2 - Optimizing RAM usage
% I noticed on the big parameter study that RAM was ticking up fast because
% we were saving the output (a,t,ct) of each simulation. We needed this data
% to plot the traces, and also to analyse the period. In version 2, we 
% now have the option to turn off these large data holds. When this option
% is toggled off, we can no longer plot the traces. See plottraces toggle.

% in order to do the data analysis without saving the data, I had to move
% the analysis into the aquisition loop. Thus the version change.

close all; clear; clc; addpath('..'); t1=tic;

plottraces = 0;

%j0array = linspace(0,0.15,60);
j0array = logspace(-3.5,-0.8239,60);
j1array = linspace(0,0.07,60);

tmax = 3000;
fig = [0 0 0 0 0 0 0 0 0];

Niter = length(j0array)*length(j1array);
[swpXj0,swpYj1] = meshgrid(j0array,j1array);
sims = struct;
Nj0 = length(j0array);
Nj1 = length(j1array);

eus_per = zeros(Nj1,Nj0);
sig_per = zeros(Nj1,Nj0);

% data aquisition (& analysis)
disp('Beginning data aquisition. AKA parameter study')

for k = 1:Niter
    disp(['Beginning simulation number ' mat2str(k) ' of ' mat2str(Niter)])
    
    jin0 = swpXj0(k);
    jin1 = swpYj1(k);
    
    [p,o]=tabakrinzelcalcium('total',tmax,...
                             'fig',fig,...
                             'makedir',0,...
                             'jin0',jin0,...
                             'jin1',jin1,...
                             'j0j1bif',1);
    
    sims.s{k}.p = p;
    
    if plottraces
        sims.s{k}.a = o.a;
        sims.s{k}.t = o.t;
        sims.s{k}.ct = o.ct;
    end
    
    % analysis
    
    % run find peaks on "a"
    [sims.s{k}.pks,sims.s{k}.locs] = ...
                 findpeaks(o.a,'minpeakheight',mean(o.a)+1*std(o.a),...
                               'minpeakprominence',0.4,...
                               'minpeakdistance',0.4/0.001);
    
    % run event classification                     
    [sims.s{k}.eupnea_i,sims.s{k}.sigh_i] = ClassifyEvents(sims.s{k}.pks,sims.s{k}.locs,o.ct);
    
    
    if ~isempty(sims.s{k}.eupnea_i)
        % eupnea mean period
        sims.s{k}.emp = mean(diff(sims.s{k}.locs(sims.s{k}.eupnea_i)*0.001));
        % check if there was only one event
        if isnan(sims.s{k}.emp); sims.s{k}.emp=0; end
    else
        sims.s{k}.emp = NaN;
    end
    
    if ~isempty(sims.s{k}.sigh_i)
        % sigh mean period
        sims.s{k}.smp = mean(diff((sims.s{k}.locs(sims.s{k}.sigh_i)*0.001)));
        % check if there was only one event
        if isnan(sims.s{k}.smp); sims.s{k}.smp=0; end
    else
        sims.s{k}.smp = NaN;
    end
    
    eus_per(k) = sims.s{k}.emp;
    sig_per(k) = sims.s{k}.smp;
    
    disp(['Concluding simulation number ' mat2str(k)])
                   
end

%% post processing
eus_freq = 1./eus_per;
sig_freq = 1./sig_per;

%%  traces figure

if plottraces
    
    % annoying stuff i have to do for subplot... to get it to match the way we
    % index into matrices. subplot's indexing runs horizontal. Matrix indexing
    % runs verticle. So I'm making a matrix with transposed indexes.
    A = zeros(Nj1,Nj0); n=0;
    for i1 = 1:Nj1
        for i2 = 1:Nj0
            n=n+1;
            A(i1,i2)=n;
        end
    end
    
    % here is the actual call to plot
    figure(1)
    for k = 1:Niter
        subplot(Nj1,Nj0,A(k))
        hold on;
        ei = sims.s{k}.eupnea_i;
        si = sims.s{k}.sigh_i;
        plot(sims.s{k}.t,sims.s{k}.a);
        if ~isempty(ei); plot(2000+sims.s{k}.locs(ei)*0.001,sims.s{k}.pks(ei),'bo'); end
        if ~isempty(si); plot(2000+sims.s{k}.locs(si)*0.001,sims.s{k}.pks(si),'ro'); end
        title(['j0 = ' mat2str(round(sims.s{k}.p.jin0,2)) '    j1 = ' mat2str(round(sims.s{k}.p.jin1,2))])
    end
end

%%  surface figure
if 1
    figure(2)
    h=surface(swpXj0,swpYj1,eus_freq); hold on;
    xlabel('j0')
    ylabel('j1')
    title('Eupnea Frequency (Hz)')
    axis([j0array(1) j0array(end) j1array(1) j1array(end)])
    set(h,'LineStyle','none')
    set(gca,'Color','k')
    colormap jet
    colorbar
    
    figure(3)
    h=surface(swpXj0,swpYj1,sig_freq); hold on;
    title('Sigh Period (Hz)')
    xlabel('j0')
    ylabel('j1')
    axis([j0array(1) j0array(end) j1array(1) j1array(end)])
    set(h,'LineStyle','none')
    set(gca,'Color','k')
    colormap jet
    colorbar
end

%% clean up
t2 = round(toc(t1));
disp(['Elapsed time is ' mat2str(t2) ' seconds, or roughly ' mat2str(round(t2/60)) ' minutes.'])

%% extra functions

function [eupnea_i,sigh_i] = ClassifyEvents(pks,locs,ct)

    % Assume all bursts are inspiratory. Unless the ct has a negative
    % slope.
    
    Nvents = length(pks);
    
    eupnea_i = [];
    sigh_i   = [];
    
    for k = 1:Nvents
        
        event_i = locs(k);
        
        if ct(event_i-1) > ct(event_i+1) % then ct is going down
            sigh_i(end+1) = k; %#ok<*AGROW>
        else
            eupnea_i(end+1) = k;
        end
        
    end
    
    % doesn't allow ONLY sighs
    if isempty(eupnea_i)
        eupnea_i = sigh_i;
        sigh_i=[];
    end
     
end