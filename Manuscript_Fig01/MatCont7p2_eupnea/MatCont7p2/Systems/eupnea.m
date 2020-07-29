function out = eupnea
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_w,par_lambdaa,par_thetaa,par_ka,par_taua,par_thetas,par_ks,par_taus,par_thetaq,par_kq,par_tauqmax,par_tauqmin,par_thetatauq,par_ktauq)
x=par_w*kmrgd(2)*kmrgd(1)-kmrgd(3);
ainf=par_lambdaa/(1+exp(4*(par_thetaa-x)/par_ka));
sinf=1/(1+exp(4*(par_thetas-kmrgd(1))/par_ks));
qinf=1/(1+exp(4*(par_thetaq-kmrgd(1))/par_kq));
tauq=(par_tauqmax-par_tauqmin)/(1+exp(4*(par_thetatauq-kmrgd(1))/par_ktauq))+par_tauqmin;
dydt=[(ainf-kmrgd(1))/par_taua;
(sinf-kmrgd(2))/par_taus;
(qinf-kmrgd(3))/tauq;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(eupnea);
y0=[0,0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_w,par_lambdaa,par_thetaa,par_ka,par_taua,par_thetas,par_ks,par_taus,par_thetaq,par_kq,par_tauqmax,par_tauqmin,par_thetatauq,par_ktauq)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_w,par_lambdaa,par_thetaa,par_ka,par_taua,par_thetas,par_ks,par_taus,par_thetaq,par_kq,par_tauqmax,par_tauqmin,par_thetatauq,par_ktauq)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_w,par_lambdaa,par_thetaa,par_ka,par_taua,par_thetas,par_ks,par_taus,par_thetaq,par_kq,par_tauqmax,par_tauqmin,par_thetatauq,par_ktauq)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_w,par_lambdaa,par_thetaa,par_ka,par_taua,par_thetas,par_ks,par_taus,par_thetaq,par_kq,par_tauqmax,par_tauqmin,par_thetatauq,par_ktauq)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_w,par_lambdaa,par_thetaa,par_ka,par_taua,par_thetas,par_ks,par_taus,par_thetaq,par_kq,par_tauqmax,par_tauqmin,par_thetatauq,par_ktauq)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_w,par_lambdaa,par_thetaa,par_ka,par_taua,par_thetas,par_ks,par_taus,par_thetaq,par_kq,par_tauqmax,par_tauqmin,par_thetatauq,par_ktauq)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_w,par_lambdaa,par_thetaa,par_ka,par_taua,par_thetas,par_ks,par_taus,par_thetaq,par_kq,par_tauqmax,par_tauqmin,par_thetatauq,par_ktauq)
