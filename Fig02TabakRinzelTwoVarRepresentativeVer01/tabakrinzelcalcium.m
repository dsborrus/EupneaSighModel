function [ param, out ] = tabakrinzelcalcium(varargin)

% scale_eupnea and scale_calcium change default parameters - handle first

param.includec=1; % calcium handling
param.includes=1; % synaptic depression
param.includetheta=1; % cellular adaptation
param.scale_eupnea=1.5;
param.scale_calcium=0.5;

if nargin>0
    disp('Non-default scale and include parameters: ')
    for i=1:nargin/2
        if strcmp(varargin{2*i-1},'includec') ...
                || strcmp(varargin{2*i-1},'includes')  ...
                || strcmp(varargin{2*i-1},'includetheta') ...
                || strcmp(varargin{2*i-1},'scale_eupnea') ...
                || strcmp(varargin{2*i-1},'scale_calcium') ...
                param.(varargin{2*i-1})=varargin{2*i};
            disp(['   ' varargin{2*i-1} ' = ' num2str(varargin{2*i}) ])
        end
    end
end

scale_eupnea = param.scale_eupnea;
scale_calcium = param.scale_calcium;

% now the default numerical parameters
param.trans = 2000;
param.total = 2200;
param.dt=0.001;
param.fig=[1 1 0 0 0 0 0 0 0 ];
param.seed=-1;
param.writetraj=1;
param.writehist=1;
param.filenametraj='default'; % default will be "run" traj (see below)
param.filenamehist='default'; % default will be "run" hist (see below)
param.thintraj=1; % 1=no thinning, 10=every 10th point, etc.
param.a_iei_binwidth=0.5;
param.ct_iei_binwidth=0.5;

% now the remainder of the default model parameters
% period should be about 4 seconds.  Active phase about 0.5 seconds.
param.N = 100; % 100
param.w=1;
param.taua=scale_eupnea*0.1; % 0.05 or 0.1
param.ka=0.2;
if param.includetheta
    param.thetaa=-0.3; % -0.3 thetaa=theta (dynamic) -0.4 to +0.1 are reasonable
    param.tausmax=scale_eupnea*0.5; % 0.5 was 0.2
    param.tausmin=scale_eupnea*0.5; % 0.5 was 0.2
else
    param.thetaa=0.18; % 0.18
    param.tausmax=scale_eupnea*5;
    param.tausmin=scale_eupnea*0.1;
end
param.thetas=0.14;
param.ks=-0.08;

param.thetataus=0.3; % 0.3
param.ktaus=-0.5; % -0.5

param.thetatheta=0.15; %  0.15+param.thetaa
param.ktheta=0.2; % 0.2

param.tauthetamax=scale_eupnea*4; % 4 is default for s&theta
param.tauthetamin=scale_eupnea*0.1; % 0.1 is default for s&theta
param.thetatautheta=0.3; % 0.3 + param.thetaa
param.ktautheta=-0.5; % -0.5

param.a_thresh=0.4; % threshold crossing for activity
param.a_ieimin=0.2; %  
param.ct_thresh=1.2; % uM; threshold crossing for total calcium
param.ct_ieimin=0;

param.v1=scale_calcium*40;
param.v2=scale_calcium*0.5;
param.v3=scale_calcium*120;
param.k3=0.3;
param.n3=2;
param.lambda=0.15;
param.thetam=0.25;
param.km=0.04;
param.thetah=0.3;
param.kh=-0.06;

param.jin0=scale_calcium*0.018; % 0.02 0.018 0
param.jin1=scale_calcium*0.04; % 0.02 0.04 0.2
param.v4=scale_calcium*0.8;
param.k4=0.3;
param.n4=4;

param.thetac=0.35; % 0.35
param.kc=0.05; % 0.05
param.lambdaa=1;
param.lambdac=1.5;

param.init_a=0.5;
param.init_s=0;
param.init_theta=0;
param.init_c=0.1;
param.init_ct=2;

param.jump_a=0;
param.jump_t=Inf;
 
if nargin>0
    disp('Non-default parameters: ')
    for i=1:nargin/2
        param.(varargin{2*i-1})=varargin{2*i};
        disp(['   ' varargin{2*i-1} ' = ' num2str(varargin{2*i}) ])
    end
end

trans=param.trans;
total=param.total;
dt=param.dt;
fig=param.fig;
seed=param.seed;

writetraj=param.writetraj;
writehist=param.writehist;
filenametraj=param.filenametraj;
filenamehist=param.filenamehist;
thintraj=param.thintraj;
a_iei_binwidth=param.a_iei_binwidth;
a_iei_binwidth=param.ct_iei_binwidth;

includec=param.includec;
includes=param.includes;
includetheta=param.includetheta;

N=param.N;
w=param.w;
taua=param.taua;
thetaa=param.thetaa;
ka=param.ka;
thetas=param.thetas;
ks=param.ks;
tausmax=param.tausmax;
tausmin=param.tausmin;
thetataus=param.thetataus;
ktaus=param.ktaus;

ktheta=param.ktheta;
thetatheta=param.thetatheta;
tauthetamax=param.tauthetamax;
tauthetamin=param.tauthetamin;
thetatautheta=param.thetatautheta;
ktautheta=param.ktautheta;

a_thresh=param.a_thresh;
a_ieimin=param.a_ieimin;
ct_thresh=param.ct_thresh;
ct_ieimin=param.ct_ieimin;

v1=param.v1;
v2=param.v2;
v3=param.v3;
k3=param.k3;
n3=param.n3;
lambda=param.lambda;
thetam=param.thetam;
km=param.km;
thetah=param.thetah;
kh=param.kh;

jin0=param.jin0;
jin1=param.jin1;
v4=param.v4;
k4=param.k4;
n4=param.n4;

thetac=param.thetac;
kc=param.kc;
lambdaa=param.lambdaa;
lambdac=param.lambdac;

init_a=param.init_a;
init_s=param.init_s;
init_theta=param.init_theta;
init_c=param.init_c;
init_ct=param.init_ct;

% for two-pulse protocol 
jump_a=param.jump_a;
jump_t=param.jump_t;

close all
if seed==-1, rng('shuffle'); else rng(seed); end

% integration of ODEs
t = [0:dt:total]; % for integration of odes; times includes transient
[ a, s, theta, c, ct] = deal(zeros(size(t)));
a(1)=init_a; s(1)=init_s; theta(1)=init_theta; c(1)=init_c; ct(1)=init_ct;
a_cross = []; ct_cross = [];
for i=2:length(t)
    
    presynpresyn=w*a(i-1);
    if includes, presynpresyn=presynpresyn*s(i-1); end
    if includetheta, presynpresyn=presynpresyn-theta(i-1); end
    
    if includec
        ainf = lambdaa*xinf(presynpresyn,thetaa,ka)+lambdac*xinf(c(i-1),thetac,kc);
    else
        ainf = xinf(presynpresyn,thetaa,ka);
    end
    
    a(i)=a(i-1)+dt*( ainf - a(i-1))/taua + ...
        randn*sqrt(dt)/N/taua*( ainf*(1-a(i-1))+(1-ainf)*a(i-1) );
    
    if includes
        sinf = xinf(a(i-1),thetas,ks);
        taus = (tausmax-tausmin)*xinf(a(i-1),thetataus,ktaus)+tausmin;
        s(i)=s(i-1)+dt*( (sinf - s(i-1))/taus );
    end
    
    if includetheta
        thetainf = xinf(a(i-1),thetatheta,ktheta);
        tautheta = (tauthetamax-tauthetamin)*xinf(a(i-1),thetatautheta,ktautheta)+tauthetamin;
        theta(i)=theta(i-1)+dt*( (thetainf - theta(i-1))/tautheta );
    end
    
    if includec
        cscale = 1;
        jpm = jin0 + jin1*a(i-1) - v4*c(i-1)^n4/(k4^n4+c(i-1)^n4);
        c(i)=c(i-1)+dt*(  cscale*( (v2+v1*finf(c(i-1),thetam,km,thetah,kh))*((ct(i-1)-c(i-1))/lambda-c(i-1))   ...
            -v3*c(i-1)^n3/(k3^n3+c(i-1)^n3) ) + jpm);
        ct(i)=ct(i-1)+dt*jpm;
    end
    
    if jump_a>0 
        for j=1:length(jump_t)
        if t(i)>=jump_t(j) && t(i-1)<jump_t(j)
        a(i)=a(i)+jump_a;
        end
        end
    end
    
    if a(i)>a_thresh && a(i-1)<=a_thresh
        a_cross = [ a_cross i]; % index (dt below)
    end
    if ct(i)<ct_thresh && ct(i-1)>=ct_thresh 
        ct_cross = [ ct_cross i]; % index (dt below)
    end
end
atrans = find(t>=trans); t=t(atrans); a=a(atrans); s=s(atrans); theta=theta(atrans); c=c(atrans); ct=ct(atrans);
out.t=t; out.a=a; out.s=s; out.theta=theta; out.c=c; out.ct=ct;
a_iei = diff(a_cross)*dt; a_iei = a_iei(find(a_iei>a_ieimin)); out.a_iei=a_iei;
ct_iei = diff(ct_cross)*dt; ct_iei = ct_iei(find(ct_iei>ct_ieimin)); out.ct_iei=ct_iei;

out.mean_a_iei=mean(a_iei); out.std_a_iei=std(a_iei); out.cv_a_iei=out.std_a_iei/out.mean_a_iei;
out.mean_ct_iei=mean(ct_iei); out.std_ct_iei=std(ct_iei); out.cv_ct_iei=out.std_ct_iei/out.mean_ct_iei;


delta_s = max(s)-min(s); out.delta_s=delta_s;
delta_theta = max(theta)-min(theta); out.delta_theta=delta_theta;
delta_a = max(a)-min(a); out.delta_a=delta_a;

% make and print figures
run = [ 'fig' datestr(now,1) '-' datestr(now,'HHMMSS') ];
theparamstr = param_str(thetaa);
titlestr = { run ; theparamstr };

if fig(1)
    figure(1);
    %cscale=c-min(c); cscale=cscale/max(cscale);
    if includec
        subplot(2,1,1);
        plot(t,s,'g',t,theta,'c',t,c,'r',t,a,'b');
        legend({'s','\theta','c','a'});
        title(titlestr)
        subplot(2,1,2); plot(t,c,'r',t,ct,'k'); legend({'c','ct'});
    else
        subplot(2,1,1);
        plot(t,s,'g',t,theta,'c',t,a,'b');
        legend({'s','\theta','a'});
        title(titlestr)
    end
    print([ run '-1.png'],'-dpng')
end

if includes
    if fig(2)
        figure(2)
        aa = 0:0.001:1;
        ss = 0:0.001:1;
        if includetheta
            fixedthetarange = [ min(theta) max(theta) ];
        else
            fixedthetarange=0;
        end
        for fixedtheta = fixedthetarange
            plot((4*(fixedtheta+thetaa)-ka*log((1-aa)./aa))./(4*w*aa),aa,'r'); hold on;
        end
        plot(xinf(aa,thetas,ks),aa,'g'); hold on;
        plot(s,a,'b'); xlabel('s'); ylabel('a');
        %axis([min(s)*0.9 max(s)*1.1 min(a)*0.9 max(a)*1.1])
        axis([0 1 min(a)*0.9 max(a)*1.1])
        title(titlestr)
        print([ run '-2.png'],'-dpng')
    end
end

if includetheta
    if fig(3)
        figure(3)
        aa = 0:0.001:1;
        for fixeds = [ min(s) max(s) ];
            plot((w*aa*fixeds+ka/4*log((1-aa)./aa)),aa,'r'); hold on;
        end
        plot(xinf(aa,thetatheta,ktheta),aa,'g'); hold on;
        plot(theta,a,'b');
        plot(theta+thetaa,a,'c');
        
        %axis([min(theta)*0.9 max(theta)*1.1 min(a)*0.9 max(a)*1.1])
        title(titlestr)
        print([ run '-3.png'],'-dpng')
    end
end

if includetheta && includes
    if fig(4)
        figure(4)
        aa = 0:0.0001:1;
        ss = 0:0.0001:1;
        for fixedtheta = [min(theta)*0.9:0.01:max(theta)*1.1]
            plot3((4*(fixedtheta+thetaa)-ka*log((1-aa)./aa))./(4*w*aa),fixedtheta*ones(size(aa)),aa,'r'); hold on;
        end
        for fixeds = [min(s)*0.9:0.01:max(s)*1.1];
            plot3(fixeds*ones(size(aa)),-thetaa+(w*aa*fixeds+ka/4*log((1-aa)./aa)),aa,'r'); hold on;
        end
        plot3(s,theta,a,'b')
        xlabel('s'); ylabel('\theta');  zlabel('a')
        
        xlabel(['s, \Delta s = ' num2str(delta_s)]);
        
        ylabel(['\theta, \Delta \theta = ' num2str(delta_theta)]);
        
        zlabel(['a, \Delta a = ' num2str(delta_a)]);
        sinf = xinf(aa,thetas,ks);
        thetainf = xinf(aa,thetatheta,ktheta);
        plot3(sinf,thetainf,aa,'k');
        %for fixeds = [min(s)*0.9:0.01:max(s)*1.1];
        for fixeds = 0:0.02:1
            plot3(fixeds*ones(size(aa)),thetainf,aa,'g'); hold on;
        end
        %for fixedtheta = [min(theta)*0.9:0.01:max(theta)*1.1]
        for fixedtheta = 0:0.02:1
            plot3(sinf,fixedtheta*ones(size(aa)),aa,'c'); hold on;
        end
        axis([ 0 1 0 1 0 1])
        %axis([min(s)*0.9 max(s)*1.1 min(theta)*0.9 max(theta)*1.1 min(a)*0.9 max(a)*1.1])
        title(titlestr)
        print([ run '-4.png'],'-dpng')
    end
    
    if fig(8)
        figure(8)
        plot(s,theta,'b'); hold on
        plot(s,theta+thetaa,'c'); hold on
        delta_theta = max(theta)-min(theta);
        xlabel('s'); ylabel(['\theta, \Delta\theta = ' num2str(delta_theta)]);
        %plot([0 1]',-thetaa*[1 1]','g--')
        sinf = xinf(aa,thetas,ks);
        thetainf = xinf(aa,thetatheta,ktheta);
        plot(sinf,thetainf,'k');
        legend({'traj','traj-\theta_a'})
        title(titlestr)
        axis([ 0 1 0 1 ])
        print([ run '-8.png'],'-dpng')
    end
end

if includec
    if fig(5)
        figure(5)
        cc = 0:0.001:0.4;
        ctct = 0:0.001:3;
        css0=fzero(@(x) jin0-v4*x^n4/(k4^n4+x^n4),0.5);
        plot(css0*ones(size(ctct)),ctct,'k'); hold on;
        css1=fzero(@(x) jin0+jin1-v4*x^n4/(k4^n4+x^n4),0.5);
        plot(css1*ones(size(ctct)),ctct,'k'); hold on;
        
        for i=1:length(cc)
            ctnull(i)=fzero(@(x) (v2+v1*finf(cc(i),thetam,km,thetah,kh))*((x-cc(i))/lambda-cc(i))-v3*cc(i)^n3/(k3^n3+cc(i)^n3),0.3);
        end
        plot(cc,ctnull,'r'); hold on;
        plot(c(:),ct(:),'b-o'); hold on; xlabel('c'); ylabel('ct'); title(titlestr)
        plot(cc,max(ctct)*cc.^n4./(k4^n4+cc.^n4),'r');
        plot(cc,max(ctct)*cc.^n5./(k5^n5+cc.^n5),'c');
        print([ run '-5.png'],'-dpng')
    end
    
    if fig(9)
        figure(9)
        cc = 0:0.01:0.4;
        presynpresyn = -1:0.01:0.2;
        [Presyn,C]  =  meshgrid(presynpresyn,cc);
        Ainf = lambdaa*xinf(Presyn,thetaa,ka)+lambdac*xinf(C,thetac,kc);
        surf(presynpresyn,cc,Ainf); hold on;
        shading flat;
        presyn=w*a.*s-theta;
        ainf = lambdaa*xinf(presyn,thetaa,ka)+lambdac*xinf(c,thetac,kc);
        plot3(presyn,c,ainf,'k');
        plot3(presyn,c,a,'b');
        xlabel('presyn'); ylabel('c'); zlabel('a'); title(titlestr)
        print([ run '-9.png'],'-dpng')
    end
    
end

if includes || includetheta
    if fig(6)
        figure(6);
        subplot(2,2,1)
        aa = 0:0.0001:1;
        tausaa = (tausmax-tausmin)*xinf(aa,thetataus,ktaus)+tausmin;
        tauthetaaa = (tauthetamax-tauthetamin)*xinf(aa,thetatautheta,ktautheta)+tauthetamin;
        plot(tausaa,aa,'r',tauthetaaa,aa,'g'); xlabel('\tau_x(a)'); ylabel('a');  title(titlestr)
        legend({'\tau_s','\tau_\theta'})
        
        subplot(2,2,2)
        apresyn = -1:0.0001:1;
        ainf = xinf(apresyn,thetaa,ka);
        plot(apresyn,ainf,'r'); xlabel('presyn'); ylabel('x_\infty');  title(titlestr)
        legend({'a_\infty'})
        
        subplot(2,2,3)
        aa = 0:0.0001:1;
        sinf = xinf(aa,thetas,ks);
        thetainf = xinf(aa,thetatheta,ktheta);
        plot(sinf,aa,'g',thetainf,aa,'b'); ylabel('a'); xlabel('x_\infty');  title(titlestr)
        legend({'s_\infty','\theta_\infty'})
        
        print([ run '-6.png'],'-dpng')
    end
end

if fig(7)
    figure(7);
    subplot(2,1,1)
    histogram(a_iei,'BinWidth',a_iei_binwidth,'BinLimits',[0,50]);
    xlabel('eupnea iei'); ylabel('#'); title(titlestr);
    
    if includec
        subplot(2,1,2)
        histogram(ct_iei,'BinWidth',ct_iei_binwidth,'BinLimits',[0,400]);
        xlabel('sigh iei'); ylabel('#');
    end
    
    print([ run '-7.png'],'-dpng')
end



dirname = [ run theparamstr ];
system([ 'mkdir ' dirname ]);

if writetraj
    % write trajectory .dat files
    atrans=find(t>=trans,1,'first');
    fid = fopen('traj.dat','w');
    for k=atrans:thintraj:length(t)
        fprintf(fid,'%6.4f  ',[t(k)-trans a(k) s(k) theta(k) c(k) ct(k) ]);
        fprintf(fid,'\n');
    end
    if ~strcmp(filenametraj,'default')
        system(['cp traj.dat ' filenametraj '.dat '  ]);
    end
    system(['mv traj.dat ' dirname ]);
    fclose(fid);
end
if writehist
    % write histogram .dat files
    [h,edges] = histcounts(a_iei,'BinWidth',a_iei_binwidth,'BinLimits',[0,50]);
    h=h/sum(h); % prob dist 
    fid = fopen('hist.dat','w');
    fprintf(fid,'%6.4f  %6.4f   \n',[(edges(1:end-1)+edges(2:end))/2 ; h]);
    fclose(fid);
    if ~strcmp(filenamehist,'default')
        system(['cp hist.dat ' filenamehist '.dat '  ]);
    end
    system(['mv hist.dat ' dirname ]);
end

% wrap tings up and collect files
save([run '.mat'])
system(['mv *mat ' dirname ]);
system(['mv *png ' dirname ]);
system(['cp tabakrinzelcalcium.m ' dirname '/tabakrinzelcalcium_script_used.txt' ]);

return

function y = xinf(x,theta,k)
y=1./(1+exp(4*(theta-x)/k));
return

function f = finf(x,theta1,k1,theta2,k2)
f=1./(1+exp((theta1-x)/k1))./(1+exp((theta2-x)/k2));
return

function str = param_str(varargin)
gap0 = '+';
str =[];
for i=1:nargin
    if i==1, gap=''; else gap=gap0; end
    param_str = inputname(i);
    param_val = varargin{i};
    str = [ str gap param_str '=' num2str(param_val) ];
end

return





