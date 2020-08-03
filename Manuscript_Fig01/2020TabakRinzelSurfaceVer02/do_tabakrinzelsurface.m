% Surface plot figure 
close all; clc; clear
system('rm -r fig*')

% These first two calls to tabakrinzelcalcium get trajectories that are 
% superimposed on the surface plot

[param, out_relax] = tabakrinzelcalcium('writetraj',1,'filenametraj','figdata1',...
    'includec',0,'includetheta',1,...
    'N',1000,...
    'scale_eupnea',1.5,...
    'trans',2000,'total',2025,'thintraj',0,'fig',[0 0 0 0 0 0 0 0 0 ],...
    'seed',0,'taua',0.001);

pause(2);

[param, out_normal] = tabakrinzelcalcium('writetraj',1,'filenametraj','figdata1',...
    'includec',0,'includetheta',1,...
    'N',1000,...
    'scale_eupnea',1.5,...
    'trans',2000,'total',2025,'thintraj',0,'fig',[0 0 0 0 0 0 0 0 0 ],...
    'seed',0);

% Make mesh and surface function 
w=1;
ka=0.2;
thetaa=-0.3;  

xinf = @(x,theta,k) 1./(1+exp(4*(theta-x)/k));
squash = @(x) 1./(1+exp(x));

a = -20:0.1:20; % 0.1
theta = 0.4:0.01:1; % 0.01
[A,THETA]=meshgrid(a,theta);
%S = (4*(THETA+thetaa)-ka*log((1-A)./A))./(4*w*A); % using squash is better
S = (4*(THETA+thetaa)-ka*A)./squash(A)/(4*w);

% The plot itself
surf(S,THETA,squash(A))
xlabel('s'); ylabel('\theta');  zlabel('a')
axis([ 0.4 1.1 -Inf Inf 0 1])
shading interp
alpha 0.7
colormap cool
hold on;

plot3(out_relax.s,out_relax.theta,out_relax.a,'-c')
plot3(out_normal.s,out_normal.theta,out_normal.a,'-b','LineWidth',1)

print('fig1.png','-dpng')
print('fig1.pdf','-dpdf')
 

% NO MORE WRITE TO DATA FILE (COULD NOT GET TIKZ TO WORK WELL)
% A1=A(:);
% THETA1=THETA(:);
% S1=S(:);
% fid = fopen('figdata2.dat','w');
% fprintf(fid,'%6.8f  %6.8f  %6.8f\n',[squash(A1) THETA1 S1]');
% fclose(fid);










