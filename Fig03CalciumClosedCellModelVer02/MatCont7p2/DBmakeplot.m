% make calcium dynamics bifurcation diagram with j0 as bif param

clear; clc;

eq = load('Systems/sigh_bif_control/diagram/Step2_findequill_hopfs.mat');
eq_filler = load('Systems\sigh_bif_control\diagram\EP_EP(2).mat');
lc = load('Systems\sigh_bif_control\diagram\Step3_LimitCycle.mat');

% equill line
j0_eq = eq.x(3,:);
css = eq.x(1,:);

hopfs = [16,39];

% equll filler line
j0_eqfill = eq_filler.x(3,:);
css_eqfill = eq_filler.x(1,:);

% limit cycle
c_lc = lc.x(1:2:end-2,:);
j0_lc = lc.x(324,:);
period = lc.x(323,:);


% figures
figure
subplot(2,1,1)
plot(j0_eq(1:hopfs(1)),css(1:hopfs(1)),'k'); hold on;
plot(j0_eq(hopfs(1):hopfs(2)),css(hopfs(1):hopfs(2)),'k--');
plot(j0_eq(hopfs(2):end),css(hopfs(2):end),'k');
plot(j0_eq(hopfs),css(hopfs),'bo')
plot(j0_eqfill,css_eqfill,'k')
xlim([0 0.2])
% idk
temp = 475;
plot(j0_lc(1:temp),max(c_lc(:,1:temp)),'r');
plot(j0_lc(1:temp),min(c_lc(:,1:temp)),'r');
ylabel('c ss')

subplot(2,1,2)
plot(j0_lc(2:temp),period(2:temp),'linewidth',2)
xlim([0 0.2])
xlabel('jin0')
ylabel('period')


if 1
    % save for LaTeX
    
end