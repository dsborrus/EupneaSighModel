% 
clear; close all; clc; myt=tic;

total = 2500;
v1_control = 20;
v1_NMB     = 300;

thetam =  0.25;   
km     =  0.04;
thetah =  0.30;
kh     = -0.06;

[p1,o1] = calciumsystem('total',total,'v1',v1_control);
[p2,o2] = calciumsystem('total',total,'v1',v1_NMB);

t1 = o1.t;   t2=o2.t;
c1 = o1.c;   c2=o2.c;
ct1= o1.ct; ct2=o2.ct;
fr1= 1/mean(o1.ct_iei); fr2= 1/mean(o2.ct_iei);

%%
% $\frac{dc}{dt} = \left[ v_1 f_\infty(c) + v_2 \right] \left[ c_{er} - c \right] - \frac{v_3 c^2}{k^2_3+c^2} + j_0 - \frac{v_4c^4}{k^4_4 + c^4}$

% figure to show what v1 f inf looks like, and how it might change
figure('position',[100 400 800 300])
subplot(1,2,1); hold on;
c = 0:0.01:1;
finf_control = v1_control*(1./(1+exp((thetam-c)/km))) .* (1./(1+exp((thetah-c)/kh)));
finf_NMB = v1_NMB*(1./(1+exp((thetam-c)/km))) .* (1./(1+exp((thetah-c)/kh)));
plot(c,finf_control,'b',c,finf_NMB,'m','linewidth',2)
legend('control','v1 doubled')
title('v1*finf(c)')
xlabel('c')
subplot(1,2,2); hold on;
fhalf1 = (1./(1+exp((thetam-c)/km)));
fhalf2 = (1./(1+exp((thetah-c)/kh)));
plot(c,fhalf1,'k',c,fhalf2,'r','linewidth',2)
title('finf(c) halves')
xlabel('c')

saveas(gcf,'finf.png')

% figure to plot trajectories
figure('position',[900 400 800 600])
subplot(2,1,1); hold on;
plot(t1,c1,'b',t1,ct1,'b--','linewidth',2)
plot(o1.locs,o1.pks,'ro')
legend('cyt Ca','total Ca')
title({['Control conditions. v1 = ' mat2str(v1_control) '.'],...
       ['sigh frequency is ' mat2str(round(fr1,2)) ' Hz.']})

subplot(2,1,2); hold on;
plot(t2,c2,'m',t2,ct2,'m--','linewidth',2)
plot(o2.locs,o2.pks,'ro')
legend('cyt Ca','total Ca')
title({['Bombesin-like conditions. v1 = ' mat2str(v1_NMB) '.'],...
       ['sigh frequency is ' mat2str(round(fr2,2)) ' Hz.']})
   
saveas(gcf,'trajectories.png')
   
toc(myt);


