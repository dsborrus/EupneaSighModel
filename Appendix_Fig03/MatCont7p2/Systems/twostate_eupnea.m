function out = twostate_eupnea
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
function dydt = fun_eval(t,kmrgd,par_w,par_thetaa,par_ka,par_taua,par_thetas,par_ks,par_taus)
ainf=1/(1+exp(4*(par_thetaa-(par_w*kmrgd(2)*kmrgd(1)))/par_ka));;
sinf=1/(1+exp(4*(par_thetas-(kmrgd(1)))/par_ks));;
dydt=[(ainf-kmrgd(1))/par_taua;;
(sinf-kmrgd(2))/par_taus;;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(twostate_eupnea);
y0=[0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_w,par_thetaa,par_ka,par_taua,par_thetas,par_ks,par_taus)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_w,par_thetaa,par_ka,par_taua,par_thetas,par_ks,par_taus)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_w,par_thetaa,par_ka,par_taua,par_thetas,par_ks,par_taus)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_w,par_thetaa,par_ka,par_taua,par_thetas,par_ks,par_taus)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_w,par_thetaa,par_ka,par_taua,par_thetas,par_ks,par_taus)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_w,par_thetaa,par_ka,par_taua,par_thetas,par_ks,par_taus)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_w,par_thetaa,par_ka,par_taua,par_thetas,par_ks,par_taus)
