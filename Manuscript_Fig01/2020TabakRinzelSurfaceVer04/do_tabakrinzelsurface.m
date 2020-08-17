% Surface plot figure 
close all; clc; clear
system('rm -r fig*')

%% Calls to simulation %%
% These first two calls to tabakrinzelcalcium get trajectories that are 
% superimposed on the surface plot

[param, out_relax] = tabakrinzelcalcium('writetraj',1,'filenametraj','figdata1',...
    'includec',0,'includetheta',1,...
    'N',inf,...
    'scale_eupnea',1.5,...
    'trans',2000,'total',2025,'thintraj',0,'fig',[0 0 0 0 0 0 0 0 0 ],...
    'seed',0,'taua',0.001);

pause(1);

[param, out_normal] = tabakrinzelcalcium('writetraj',1,'filenametraj','figdata1',...
    'includec',0,'includetheta',1,...
    'N',inf,...
    'scale_eupnea',1.5,...
    'trans',2000,'total',2025,'thintraj',0,'fig',[0 0 0 0 0 0 0 0 0 ],...
    'seed',0);

%% Make mesh and surface function %%
w=1;
ka=2; %0.2
thetaa=0; %-0.3
lambdaa=10;
amin = 1;

squash = @(x) 1./(1+exp(x));
squash2 = @(x) (10+exp(x)*1)./(1+exp(x));

a = -20:0.01:20; % 0.1
theta = 1:0.1:7.5; % 0.01
[A,THETA]=meshgrid(a,theta);
%S = (4*(THETA+thetaa)-ka*log((1-A)./A))./(4*w*A); % using squash is better
%S = (4*(THETA+thetaa)-ka*A)./squash(A)/(4*w);
%S = abs( (-(ka/4)*log((lambdaa-A)./(A-amin))+THETA+thetaa)./(w*A));
%S = abs( (-(ka/4)*log((lambdaa-A)./(A-amin))+THETA+thetaa)./(w*A));
S = (4*(THETA+thetaa)-ka*A)./squash2(A)/(4*w);

%% The activity plot itself %%
f=figure;%figure('position',[400 400 1100 800])
surf(S,THETA,squash2(A))
%xlabel('s'); ylabel('\theta');  zlabel('a')
axis([ 0.4 1.0 1 7.5 0 lambdaa])%axis([ 0.4 1.1 1 7.5 0 lambdaa])
view(-49,17)
%view(-51,15)
shading interp
alpha 0.7
colormap parula
hold on;
xticks([0.4 1]); yticks([1 7.5]); zticks([0 10])
set(gca,'fontsize',18)

p1=2557;
p2=2615;
p3=2690;
p4=4000;

ms=10;

plot3(out_relax.s,out_relax.theta,out_relax.a,'-c','LineWidth',2);
% placing directional arrows on the trajectory. This effect only works from
% the default viewing angle!!!***
plot3(out_relax.s(p1),out_relax.theta(p1),out_relax.a(p1),'c','marker','^','markerfacecolor','c','markersize',ms,'HandleVisibility','off')
plot3(out_relax.s(p2),out_relax.theta(p2),out_relax.a(p2),'c','marker','<','markerfacecolor','c','markersize',ms,'HandleVisibility','off')
plot3(out_relax.s(p3),out_relax.theta(p3),out_relax.a(p3),'c','marker','v','markerfacecolor','c','markersize',ms,'HandleVisibility','off')
plot3(out_relax.s(p4),out_relax.theta(p4),out_relax.a(p4),'c','marker','>','markerfacecolor','c','markersize',ms,'HandleVisibility','off')

p1=600;
p2=1250;
p3=5150;
p4=5400;

plot3(out_normal.s,out_normal.theta,out_normal.a,'-b','LineWidth',2);
% placing directional arrows on the trajectory. This effect only works from
% the default viewing angle!!!***
%plot3(out_normal.s(p1),out_normal.theta(p1),out_normal.a(p1),'b','marker','>','markerfacecolor','b','markersize',ms,'HandleVisibility','off')
plot3(out_normal.s(p2),out_normal.theta(p2),out_normal.a(p2),'b','marker','>','markerfacecolor','b','markersize',ms,'HandleVisibility','off')
plot3(out_normal.s(p3),out_normal.theta(p3),out_normal.a(p3),'b','marker','>','markerfacecolor','b','markersize',ms,'HandleVisibility','off')
plot3(out_normal.s(p4),out_normal.theta(p4),out_normal.a(p4),'b','marker','>','markerfacecolor','b','markersize',ms,'HandleVisibility','off')

%% Background nulc line traces %%

aa = linspace(0,lambdaa,1e3);
% s nulcline
thetas = param.thetas;
ks = param.ks;
sinf = xinf(aa,thetas,ks);
temptheta = zeros(1e3,1)+max(max(THETA));
plot3(sinf,temptheta,aa,'color',[0.7 0 0.7],'MarkerIndices',1:20:800);

% theta nulcline
thetatheta = param.thetatheta;
ktheta = param.ktheta;
thetainf = lambdaa*xinf(aa,thetatheta,ktheta);
temps = zeros(1e3,1)+f.CurrentAxes.XLim(2);
plot3(temps,thetainf,aa,'color',[0 0.5 0],'marker','none','MarkerIndices',1:10:length(temps))

% legend
l=legend('','','','','');
axP = get(gca,'Position');
l.Position=[0.7979 0.3602 0.1254 0.4048];
legend boxoff
set(gca, 'Position', axP);


%% Output %%
print('fig1.png','-dpng')
print('fig1.pdf','-dpdf')
print('../tikz/PanelB/PanelB.png','-dpng')
 

% NO MORE WRITE TO DATA FILE (COULD NOT GET TIKZ TO WORK WELL)
% A1=A(:);
% THETA1=THETA(:);
% S1=S(:);
% fid = fopen('figdata2.dat','w');
% fprintf(fid,'%6.8f  %6.8f  %6.8f\n',[squash(A1) THETA1 S1]');
% fclose(fid);




function y = xinf(x,theta,k)
y=1./(1+exp(4*(theta-x)/k));
end



