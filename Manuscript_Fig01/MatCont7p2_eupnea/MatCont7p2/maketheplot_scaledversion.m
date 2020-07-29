clear; close all; clc;


hps = load('Systems/ca_sys_ver2(scaled)/diagram/Hopf_bif_ver2.mat');

% hopf figure
v1 = hps.x(3,:);
c  = hps.x(1,:);
poi1 = hps.s(2).index; poi2 = hps.s(3).index;
hb1 = [hps.x(3,poi1) hps.x(1,poi1)];
hb2 = [hps.x(3,poi2) hps.x(1,poi2)];

lc = load('Systems/ca_sys_ver2(scaled)/diagram/LimitCycle_ver2.mat');

% period figure
v1p = lc.x(end,:);
p  = lc.x(end-1,:);
f  = 1./p;

% limit cycle
LC = (lc.x(1:2:end-2,:));

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
    
    plot(v1p,max(LC),'b','linewidth',2)
    plot(v1p,min(LC),'b','linewidth',2,'HandleVisibility','off')
    set(gca,'fontsize',11)
    
    legend('sss of ca','uss of ca','hopf bifurcation','limit cycle upper and lower bounds')
    
    subplot(2,1,2); hold on;
    plot(v1p,f,'r','linewidth',2)
    xlim([min(v1) max(v1)])
    xlabel('v1');
    ylabel('freq (Hz)');
    title('frequency of oscillations (Hz)')
    set(gca,'fontsize',11)
end

biLC_v1p = v1p;
biLC_max = max(LC);
biLC_min = min(LC);

save('bi_LC.mat','biLC_v1p', 'biLC_max', 'biLC_min')
