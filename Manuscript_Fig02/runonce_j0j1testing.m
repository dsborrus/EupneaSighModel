% this script just runs the model once with speccficied j0 and j1 values
clear; close all

tau_ct = 2.5;

j0 = 0.005;
j1 = 0.00512; %0.04;
v4 = 0.16;
tauthetamax = 8;% default is 6

ct_thresh = 1.3;

scale_calcium = 0.5;
amax = 10;

%jin0=scale_calcium      *   j0 / tau_ct; % 0.02 0.018 0
%jin1=scale_calcium/amax *   j1 / tau_ct; % 0.02 0.04 0.2
%v4_ =scale_calcium      *   v4 / tau_ct;

jin0 = j0;
jin1 = j1;
v4_ = v4;
[p,out] = tabakrinzelcalcium('jin1',jin1,'jin0',jin0,...
                             'v4',v4_,...
                             'thetaa',-2,...
                             'tauthetamax',tauthetamax,...
                             'fig',[0 0 0 0 1 0 0 0 0 0],...
                             'total',3000,...
                             'ct_thresh',ct_thresh);

figure
ax1=subplot(2,1,1); hold on;
plot(out.t,out.a,out.t,out.s,out.t,out.theta)       

ax2=subplot(2,1,2); hold on;
plot(out.t,out.c,out.t,out.ct)
legend('c','ct')

linkaxes([ax1 ax2],'x')

disp(['sigh period avg is ' mat2str(out.mean_ct_iei) ' sec'])