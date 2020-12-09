%% Manuscript Figure #03 - NMB Experiments
% script to make and export figures for figure 3 of manuscript
clear; close all; clc;

%% Simulation
myt=tic;

total = 2500;
v1_control = 20;
v1_NMB     = 40;

thetam =  0.25;   
km     =  0.04;
thetah =  0.30;
kh     = -0.06;

[p1,o1] = tabakrinzelcalcium('total',total,'v1',v1_control);
[p2,o2] = tabakrinzelcalcium('total',total,'v1',v1_NMB);

t1 = o1.t;   t2=o2.t;
c1 = o1.c;   c2=o2.c;
ct1= o1.ct; ct2=o2.ct;
fr1= 1/mean(o1.ct_iei); fr2= 1/mean(o2.ct_iei);

                                 
%% Summary plot

% figure to show what v1 f inf looks like, and how it might change
figure('position',[100 400 800 300])
subplot(1,2,1); hold on;
c = 0:0.01:1;
finf_control = v1_control*(1./(1+exp((thetam-c)/km))) .* (1./(1+exp((thetah-c)/kh)));
finf_NMB = v1_NMB*(1./(1+exp((thetam-c)/km))) .* (1./(1+exp((thetah-c)/kh)));
plot(c,finf_control,'b',c,finf_NMB,'m','linewidth',2)
legend('control (v1=20)','v1 doubled (v1=40)')
title('v1*finf(c)')
xlabel('c')
subplot(1,2,2); hold on;
fhalf1 = (1./(1+exp((thetam-c)/km)));
fhalf2 = (1./(1+exp((thetah-c)/kh)));
plot(c,fhalf1,'k',c,fhalf2,'r','linewidth',2)
title('finf(c) halves')
xlabel('c')

% figure to plot trajectories
figure('position',[900 400 800 600])
subplot(2,1,1); hold on;
plot(t1,c1,'b',t1,ct1,'b--','linewidth',2)
%plot(o1.locs,o1.pks,'ro')
legend('cyt Ca','total Ca')
title({['Control conditions. v1 = ' mat2str(v1_control) '.'],...
       ['sigh frequency is ' mat2str(round(fr1,2)) ' Hz.']})

subplot(2,1,2); hold on;
plot(t2,c2,'m',t2,ct2,'m--','linewidth',2)
%plot(o2.locs,o2.pks,'ro')
legend('cyt Ca','total Ca')
title({['Bombesin-like conditions. v1 = ' mat2str(v1_NMB) '.'],...
       ['sigh frequency is ' mat2str(round(fr2,2)) ' Hz.']})


% figure to plot nulc lines

figure
cc = 0:0.001:0.4;
v2=0.25;v3=60;k3=0.3;n3=2;lambda=0.15;thetam=0.25;km=0.04;thetah=0.3;kh=-0.06;j0=0.009;v4=0.4;k4=0.3;n4=4;

for i=1:length(cc)
    %ctnull1(i)=fzero(@(x) (v2+v1_control*finf(cc(i),thetam,km,thetah,kh))*((x-cc(i))/lambda-cc(i))-v3*cc(i)^n3/(k3^n3+cc(i)^n3),0.3);
    ctnull_con(i)=fzero(@(x) (v2+v1_control*finf(cc(i),thetam,km,thetah,kh))*((x-cc(i))/lambda-cc(i))-v3*cc(i)^n3/(k3^n3+cc(i)^n3)+j0-((v4*cc(i)^4)/(k4+cc(i)^4)),0.3);
    ctnull_nmb(i)=fzero(@(x) (v2+v1_NMB*finf(cc(i),thetam,km,thetah,kh))*((x-cc(i))/lambda-cc(i))-v3*cc(i)^n3/(k3^n3+cc(i)^n3)+j0-((v4*cc(i)^4)/(k4+cc(i)^4)),0.3);    
    ctnull_10(i)=fzero(@(x) (v2+10*finf(cc(i),thetam,km,thetah,kh))*((x-cc(i))/lambda-cc(i))-v3*cc(i)^n3/(k3^n3+cc(i)^n3)+j0-((v4*cc(i)^4)/(k4+cc(i)^4)),0.3);    
    %ctnull_80(i)=fzero(@(x) (v2+80*finf(cc(i),thetam,km,thetah,kh))*((x-cc(i))/lambda-cc(i))-v3*cc(i)^n3/(k3^n3+cc(i)^n3)+j0-((v4*cc(i)^4)/(k4+cc(i)^4)),0.3);
end

%plot(ctnull1,'b'); hold on; plot(ctnull2,'r')

plot(cc,ctnull_10,'color',[0.4 0 0.4],'linewidth',2); hold on;
plot(cc,ctnull_con,'r','linewidth',2);
plot(cc,ctnull_nmb,'m','linewidth',2);
%plot(cc,ctnull_80,'color',[0.9 0.1 0.9],'linewidth',2);
plot( [((j0*k4^4)/(v4-j0))^(1/4) ((j0*k4^4)/(v4-j0))^(1/4)], [0 max([ctnull_10 ctnull_con ctnull_nmb])],'b','linewidth',2 )      

legend("c'=0 (v1=10)","c'=0 control (v1=20)","c'=0 nmb (v1=40)","ct'=0 nullcline")
xlabel('c')
ylabel('ctot')
title('Phase diagram, Ca subsystem')
axis tight

%% Export data
res = 50;
control_traj = [t1(1:res:end)-2000; o1.a(1:res:end); o1.c(1:res:end); o1.ct(1:res:end)]';
save('tikz/data/control_traj.dat','control_traj','-ascii');
NMB_traj = [t2(1:res:end)-2000; o2.a(1:res:end); o2.c(1:res:end); o2.ct(1:res:end)]';
save('tikz/data/NMB_traj.dat','NMB_traj','-ascii');

% save parameters to tex file

fileID = fopen('tikz/PanelA/defs.tex' ,'w');
% number of sighs
fprintf(fileID,['\\def \\maxAcontrol{' num2str(max(o1.a)) '}\n']);
fprintf(fileID,['\\def \\maxAnmb{' num2str(max(o2.a)) '}\n']);
fprintf(fileID,['\\def \\cmin{' num2str(min(o1.c)) '}\n']);
fprintf(fileID,['\\def \\cmax{' num2str(max(o1.c)) '}\n']);
fprintf(fileID,['\\def \\ctmin{' num2str(min(o1.ct)) '}\n']);
fprintf(fileID,['\\def \\ctmax{' num2str(max(o1.ct)) '}\n']);
fprintf(fileID,['\\def \\ctminNMB{' num2str(min(o2.ct)) '}\n']);
fprintf(fileID,['\\def \\ctmaxNMB{' num2str(max(o2.ct)) '}\n']);


%% Additional Functions

function y = xinf(x,theta,k,lambda)
    y=lambda./(1+exp(4*(theta-x)/k));
end

function f = finf(x,theta1,k1,theta2,k2)
    f=1./(1+exp((theta1-x)/k1))./(1+exp((theta2-x)/k2));
end
