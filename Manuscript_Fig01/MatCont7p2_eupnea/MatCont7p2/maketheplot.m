clear; close all; clc;


hps = load('Systems/ca_sys/diagram/EP_EP(2)_HopfEquillCurve.mat');

% hopf figure
v1 = hps.x(3,:);
c  = hps.x(1,:);
poi1 = hps.s(2).index; poi2 = hps.s(3).index;
hb1 = [hps.x(3,poi1) hps.x(1,poi1)];
hb2 = [hps.x(3,poi2) hps.x(1,poi2)];

lc1 = load('Systems/ca_sys/diagram/limit_cycle_startmaybe.mat');
lc2 = load('Systems/ca_sys/diagram/limit_cycle_finalmaybe.mat');

% period figure
v11 = lc1.x(end,:);   v12 = fliplr(lc2.x(end,:));
p1  = lc1.x(end-1,:); p2  = fliplr(lc2.x(end-1,:));
f1  = 1./p1;           f2  = 1./p2;

% limit cycle
LC1 = fliplr(lc1.x(1:2:end-2,:));
LC2 = fliplr(lc2.x(1:2:end-2,:));

if 1
    figure;
    
    subplot(2,1,1); hold on; 
    title('Ca oscillations as a function of maximum Ca ER efflux (v1)')
    
    plot(v1(1:poi1),c(1:poi1),'k','linewidth',2); 
    plot(v1(poi1:poi2),c(poi1:poi2),'k--','linewidth',2);
    plot(v1(poi2:end),c(poi2:end),'k','linewidth',2,'HandleVisibility','off');
    
    plot(hb1(1),hb1(2),'ko','markersize',11,'linewidth',2)
    plot(hb2(1),hb2(2),'ko','markersize',11,'linewidth',2,'HandleVisibility','off')
    xlim([min(v1) max(v1)])
    ylim([0 0.6])
    xlabel('v1'); ylabel('c');
    
    plot(v11,max(LC1),'b','linewidth',2)
    plot(v11,min(LC1),'b','linewidth',2,'HandleVisibility','off')
    plot(v12,max(LC2),'b','linewidth',2,'HandleVisibility','off')
    plot(v12,min(LC2),'b','linewidth',2,'HandleVisibility','off')
    set(gca,'fontsize',11)
    
    legend('sss of ca','uss of ca','hopf bifurcation','limit cycle upper and lower bounds')
    
    
    subplot(2,1,2); hold on;
    plot(v11,f1,'r','linewidth',2)
    plot(v12,f2,'r','linewidth',2)
    xlim([min(v1) max(v1)])
    xlabel('v1');
    ylabel('freq (Hz)');
    title('frequency of oscillations (Hz)')
    set(gca,'fontsize',11)
end