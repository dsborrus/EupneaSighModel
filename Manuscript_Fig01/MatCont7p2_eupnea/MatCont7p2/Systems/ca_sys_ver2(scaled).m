function out = ca_sys_ver2(scaled)
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
function dydt = fun_eval(t,kmrgd,par_v1,par_v2,par_v3,par_k3,par_n3,par_lambda,par_thetam,par_km,par_thetah,par_kh,par_jin0,par_v4,par_k4,par_n4,par_eta)
finf=(1/(1+exp((par_thetam-kmrgd(1))/par_km)))*(1/(1+exp((par_thetah-kmrgd(1))/par_kh)));;
exfl=par_jin0-(par_v4*kmrgd(1)^4)/(par_k4^4+kmrgd(1)^4);;
cer=(kmrgd(2)-kmrgd(1))/par_lambda;;
dydt=[par_eta*((par_v1*finf+par_v2)*(cer-kmrgd(1))-(par_v3*kmrgd(1)^2)/(par_k3^2+kmrgd(1)^2))+exfl;;
exfl;;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(ca_sys_ver2(scaled));
y0=[0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_v1,par_v2,par_v3,par_k3,par_n3,par_lambda,par_thetam,par_km,par_thetah,par_kh,par_jin0,par_v4,par_k4,par_n4,par_eta)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_v1,par_v2,par_v3,par_k3,par_n3,par_lambda,par_thetam,par_km,par_thetah,par_kh,par_jin0,par_v4,par_k4,par_n4,par_eta)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_v1,par_v2,par_v3,par_k3,par_n3,par_lambda,par_thetam,par_km,par_thetah,par_kh,par_jin0,par_v4,par_k4,par_n4,par_eta)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_v1,par_v2,par_v3,par_k3,par_n3,par_lambda,par_thetam,par_km,par_thetah,par_kh,par_jin0,par_v4,par_k4,par_n4,par_eta)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_v1,par_v2,par_v3,par_k3,par_n3,par_lambda,par_thetam,par_km,par_thetah,par_kh,par_jin0,par_v4,par_k4,par_n4,par_eta)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_v1,par_v2,par_v3,par_k3,par_n3,par_lambda,par_thetam,par_km,par_thetah,par_kh,par_jin0,par_v4,par_k4,par_n4,par_eta)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_v1,par_v2,par_v3,par_k3,par_n3,par_lambda,par_thetam,par_km,par_thetah,par_kh,par_jin0,par_v4,par_k4,par_n4,par_eta)
