function out = sigh_hopf
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_v1,par_v2,par_v3,par_v4,par_k3,par_k4,par_thetam,par_km,par_thetah,par_kh,par_lambda,par_j0,par_j1,par_a)
finf=(1/(1+exp((par_thetam-kmrgd(1))/par_km)))+(1/(1+exp((par_thetah-kmrgd(1))/par_kh)));
cer=(kmrgd(2)-kmrgd(1))/par_lambda;
dydt=[(par_v1+par_v2*finf)*(cer-kmrgd(1))-(par_v3*kmrgd(1)^2)/(par_k3^2+kmrgd(1)^2)+par_j0+par_j1*par_a-(par_v4*kmrgd(1)^4)/(par_k4^4+kmrgd(1)^4);
par_j0+par_j1*par_a-(par_v4*kmrgd(1)^4)/(par_k4^4+kmrgd(1)^4);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(sigh_hopf);
y0=[0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_v1,par_v2,par_v3,par_v4,par_k3,par_k4,par_thetam,par_km,par_thetah,par_kh,par_lambda,par_j0,par_j1,par_a)
jac=[ (2*kmrgd(1)^3*par_v3)/(kmrgd(1)^2 + par_k3^2)^2 - par_v2*(kmrgd(1) + (kmrgd(1) - kmrgd(2))/par_lambda)*(exp(-(kmrgd(1) - par_thetah)/par_kh)/(par_kh*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) + exp(-(kmrgd(1) - par_thetam)/par_km)/(par_km*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2)) - (2*kmrgd(1)*par_v3)/(kmrgd(1)^2 + par_k3^2) - (1/par_lambda + 1)*(par_v1 + par_v2*(1/(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1) + 1/(exp(-(kmrgd(1) - par_thetam)/par_km) + 1))) - (4*kmrgd(1)^3*par_v4)/(kmrgd(1)^4 + par_k4^4) + (4*kmrgd(1)^7*par_v4)/(kmrgd(1)^4 + par_k4^4)^2 , (par_v1 + par_v2*(1/(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1) + 1/(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)))/par_lambda ; (4*kmrgd(1)^7*par_v4)/(kmrgd(1)^4 + par_k4^4)^2 - (4*kmrgd(1)^3*par_v4)/(kmrgd(1)^4 + par_k4^4) , 0 ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_v1,par_v2,par_v3,par_v4,par_k3,par_k4,par_thetam,par_km,par_thetah,par_kh,par_lambda,par_j0,par_j1,par_a)
jacp=[ - kmrgd(1) - (kmrgd(1) - kmrgd(2))/par_lambda , -(1/(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1) + 1/(exp(-(kmrgd(1) - par_thetam)/par_km) + 1))*(kmrgd(1) + (kmrgd(1) - kmrgd(2))/par_lambda) , -kmrgd(1)^2/(kmrgd(1)^2 + par_k3^2) , -kmrgd(1)^4/(kmrgd(1)^4 + par_k4^4) , (2*kmrgd(1)^2*par_k3*par_v3)/(kmrgd(1)^2 + par_k3^2)^2 , (4*kmrgd(1)^4*par_k4^3*par_v4)/(kmrgd(1)^4 + par_k4^4)^2 , (par_v2*exp(-(kmrgd(1) - par_thetam)/par_km)*(kmrgd(1) + (kmrgd(1) - kmrgd(2))/par_lambda))/(par_km*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2) , (par_v2*exp(-(kmrgd(1) - par_thetam)/par_km)*(kmrgd(1) + (kmrgd(1) - kmrgd(2))/par_lambda)*(kmrgd(1) - par_thetam))/(par_km^2*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2) , (par_v2*exp(-(kmrgd(1) - par_thetah)/par_kh)*(kmrgd(1) + (kmrgd(1) - kmrgd(2))/par_lambda))/(par_kh*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) , (par_v2*exp(-(kmrgd(1) - par_thetah)/par_kh)*(kmrgd(1) + (kmrgd(1) - kmrgd(2))/par_lambda)*(kmrgd(1) - par_thetah))/(par_kh^2*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) , ((kmrgd(1) - kmrgd(2))*(par_v1 + par_v2*(1/(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1) + 1/(exp(-(kmrgd(1) - par_thetam)/par_km) + 1))))/par_lambda^2 , 1 , par_a , par_j1 ; 0 , 0 , 0 , -kmrgd(1)^4/(kmrgd(1)^4 + par_k4^4) , 0 , (4*kmrgd(1)^4*par_k4^3*par_v4)/(kmrgd(1)^4 + par_k4^4)^2 , 0 , 0 , 0 , 0 , 0 , 1 , par_a , par_j1 ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_v1,par_v2,par_v3,par_v4,par_k3,par_k4,par_thetam,par_km,par_thetah,par_kh,par_lambda,par_j0,par_j1,par_a)
hess1=[ par_v2*(kmrgd(1) + (kmrgd(1) - kmrgd(2))/par_lambda)*(exp(-(kmrgd(1) - par_thetah)/par_kh)/(par_kh^2*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) - (2*exp(-(2*(kmrgd(1) - par_thetah))/par_kh))/(par_kh^2*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^3) + exp(-(kmrgd(1) - par_thetam)/par_km)/(par_km^2*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2) - (2*exp(-(2*(kmrgd(1) - par_thetam))/par_km))/(par_km^2*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^3)) - (2*par_v3)/(kmrgd(1)^2 + par_k3^2) - 2*par_v2*(1/par_lambda + 1)*(exp(-(kmrgd(1) - par_thetah)/par_kh)/(par_kh*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) + exp(-(kmrgd(1) - par_thetam)/par_km)/(par_km*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2)) + (10*kmrgd(1)^2*par_v3)/(kmrgd(1)^2 + par_k3^2)^2 - (8*kmrgd(1)^4*par_v3)/(kmrgd(1)^2 + par_k3^2)^3 - (12*kmrgd(1)^2*par_v4)/(kmrgd(1)^4 + par_k4^4) + (44*kmrgd(1)^6*par_v4)/(kmrgd(1)^4 + par_k4^4)^2 - (32*kmrgd(1)^10*par_v4)/(kmrgd(1)^4 + par_k4^4)^3 , (par_v2*(exp(-(kmrgd(1) - par_thetah)/par_kh)/(par_kh*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) + exp(-(kmrgd(1) - par_thetam)/par_km)/(par_km*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2)))/par_lambda ; (44*kmrgd(1)^6*par_v4)/(kmrgd(1)^4 + par_k4^4)^2 - (12*kmrgd(1)^2*par_v4)/(kmrgd(1)^4 + par_k4^4) - (32*kmrgd(1)^10*par_v4)/(kmrgd(1)^4 + par_k4^4)^3 , 0 ];
hess2=[ (par_v2*(exp(-(kmrgd(1) - par_thetah)/par_kh)/(par_kh*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) + exp(-(kmrgd(1) - par_thetam)/par_km)/(par_km*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2)))/par_lambda , 0 ; 0 , 0 ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_v1,par_v2,par_v3,par_v4,par_k3,par_k4,par_thetam,par_km,par_thetah,par_kh,par_lambda,par_j0,par_j1,par_a)
hessp1=[ - 1/par_lambda - 1 , 1/par_lambda ; 0 , 0 ];
hessp2=[ - (kmrgd(1) + (kmrgd(1) - kmrgd(2))/par_lambda)*(exp(-(kmrgd(1) - par_thetah)/par_kh)/(par_kh*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) + exp(-(kmrgd(1) - par_thetam)/par_km)/(par_km*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2)) - (1/(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1) + 1/(exp(-(kmrgd(1) - par_thetam)/par_km) + 1))*(1/par_lambda + 1) , (1/(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1) + 1/(exp(-(kmrgd(1) - par_thetam)/par_km) + 1))/par_lambda ; 0 , 0 ];
hessp3=[ (2*kmrgd(1)^3)/(kmrgd(1)^2 + par_k3^2)^2 - (2*kmrgd(1))/(kmrgd(1)^2 + par_k3^2) , 0 ; 0 , 0 ];
hessp4=[ (4*kmrgd(1)^7)/(kmrgd(1)^4 + par_k4^4)^2 - (4*kmrgd(1)^3)/(kmrgd(1)^4 + par_k4^4) , 0 ; (4*kmrgd(1)^7)/(kmrgd(1)^4 + par_k4^4)^2 - (4*kmrgd(1)^3)/(kmrgd(1)^4 + par_k4^4) , 0 ];
hessp5=[ (4*kmrgd(1)*par_k3*par_v3)/(kmrgd(1)^2 + par_k3^2)^2 - (8*kmrgd(1)^3*par_k3*par_v3)/(kmrgd(1)^2 + par_k3^2)^3 , 0 ; 0 , 0 ];
hessp6=[ (16*kmrgd(1)^3*par_k4^3*par_v4)/(kmrgd(1)^4 + par_k4^4)^2 - (32*kmrgd(1)^7*par_k4^3*par_v4)/(kmrgd(1)^4 + par_k4^4)^3 , 0 ; (16*kmrgd(1)^3*par_k4^3*par_v4)/(kmrgd(1)^4 + par_k4^4)^2 - (32*kmrgd(1)^7*par_k4^3*par_v4)/(kmrgd(1)^4 + par_k4^4)^3 , 0 ];
hessp7=[ (par_v2*exp(-(kmrgd(1) - par_thetam)/par_km)*(1/par_lambda + 1))/(par_km*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2) - par_v2*(kmrgd(1) + (kmrgd(1) - kmrgd(2))/par_lambda)*(exp(-(kmrgd(1) - par_thetam)/par_km)/(par_km^2*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2) - (2*exp(-(2*(kmrgd(1) - par_thetam))/par_km))/(par_km^2*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^3)) , -(par_v2*exp(-(kmrgd(1) - par_thetam)/par_km))/(par_km*par_lambda*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2) ; 0 , 0 ];
hessp8=[ par_v2*(kmrgd(1) + (kmrgd(1) - kmrgd(2))/par_lambda)*(exp(-(kmrgd(1) - par_thetam)/par_km)/(par_km^2*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2) - (exp(-(kmrgd(1) - par_thetam)/par_km)*(kmrgd(1) - par_thetam))/(par_km^3*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2) + (2*exp(-(2*(kmrgd(1) - par_thetam))/par_km)*(kmrgd(1) - par_thetam))/(par_km^3*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^3)) + (par_v2*exp(-(kmrgd(1) - par_thetam)/par_km)*(1/par_lambda + 1)*(kmrgd(1) - par_thetam))/(par_km^2*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2) , -(par_v2*exp(-(kmrgd(1) - par_thetam)/par_km)*(kmrgd(1) - par_thetam))/(par_km^2*par_lambda*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2) ; 0 , 0 ];
hessp9=[ (par_v2*exp(-(kmrgd(1) - par_thetah)/par_kh)*(1/par_lambda + 1))/(par_kh*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) - par_v2*(kmrgd(1) + (kmrgd(1) - kmrgd(2))/par_lambda)*(exp(-(kmrgd(1) - par_thetah)/par_kh)/(par_kh^2*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) - (2*exp(-(2*(kmrgd(1) - par_thetah))/par_kh))/(par_kh^2*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^3)) , -(par_v2*exp(-(kmrgd(1) - par_thetah)/par_kh))/(par_kh*par_lambda*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) ; 0 , 0 ];
hessp10=[ par_v2*(kmrgd(1) + (kmrgd(1) - kmrgd(2))/par_lambda)*(exp(-(kmrgd(1) - par_thetah)/par_kh)/(par_kh^2*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) - (exp(-(kmrgd(1) - par_thetah)/par_kh)*(kmrgd(1) - par_thetah))/(par_kh^3*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) + (2*exp(-(2*(kmrgd(1) - par_thetah))/par_kh)*(kmrgd(1) - par_thetah))/(par_kh^3*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^3)) + (par_v2*exp(-(kmrgd(1) - par_thetah)/par_kh)*(1/par_lambda + 1)*(kmrgd(1) - par_thetah))/(par_kh^2*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) , -(par_v2*exp(-(kmrgd(1) - par_thetah)/par_kh)*(kmrgd(1) - par_thetah))/(par_kh^2*par_lambda*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) ; 0 , 0 ];
hessp11=[ (par_v1 + par_v2*(1/(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1) + 1/(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)))/par_lambda^2 + (par_v2*(kmrgd(1) - kmrgd(2))*(exp(-(kmrgd(1) - par_thetah)/par_kh)/(par_kh*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) + exp(-(kmrgd(1) - par_thetam)/par_km)/(par_km*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2)))/par_lambda^2 , -(par_v1 + par_v2*(1/(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1) + 1/(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)))/par_lambda^2 ; 0 , 0 ];
hessp12=[ 0 , 0 ; 0 , 0 ];
hessp13=[ 0 , 0 ; 0 , 0 ];
hessp14=[ 0 , 0 ; 0 , 0 ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
hessp(:,:,5) =hessp5;
hessp(:,:,6) =hessp6;
hessp(:,:,7) =hessp7;
hessp(:,:,8) =hessp8;
hessp(:,:,9) =hessp9;
hessp(:,:,10) =hessp10;
hessp(:,:,11) =hessp11;
hessp(:,:,12) =hessp12;
hessp(:,:,13) =hessp13;
hessp(:,:,14) =hessp14;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_v1,par_v2,par_v3,par_v4,par_k3,par_k4,par_thetam,par_km,par_thetah,par_kh,par_lambda,par_j0,par_j1,par_a)
tens31=[ 3*par_v2*(1/par_lambda + 1)*(exp(-(kmrgd(1) - par_thetah)/par_kh)/(par_kh^2*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) - (2*exp(-(2*(kmrgd(1) - par_thetah))/par_kh))/(par_kh^2*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^3) + exp(-(kmrgd(1) - par_thetam)/par_km)/(par_km^2*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2) - (2*exp(-(2*(kmrgd(1) - par_thetam))/par_km))/(par_km^2*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^3)) - par_v2*(kmrgd(1) + (kmrgd(1) - kmrgd(2))/par_lambda)*(exp(-(kmrgd(1) - par_thetah)/par_kh)/(par_kh^3*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) - (6*exp(-(2*(kmrgd(1) - par_thetah))/par_kh))/(par_kh^3*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^3) + (6*exp(-(3*(kmrgd(1) - par_thetah))/par_kh))/(par_kh^3*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^4) + exp(-(kmrgd(1) - par_thetam)/par_km)/(par_km^3*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2) - (6*exp(-(2*(kmrgd(1) - par_thetam))/par_km))/(par_km^3*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^3) + (6*exp(-(3*(kmrgd(1) - par_thetam))/par_km))/(par_km^3*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^4)) + (24*kmrgd(1)*par_v3)/(kmrgd(1)^2 + par_k3^2)^2 - (24*kmrgd(1)*par_v4)/(kmrgd(1)^4 + par_k4^4) - (72*kmrgd(1)^3*par_v3)/(kmrgd(1)^2 + par_k3^2)^3 + (48*kmrgd(1)^5*par_v3)/(kmrgd(1)^2 + par_k3^2)^4 + (312*kmrgd(1)^5*par_v4)/(kmrgd(1)^4 + par_k4^4)^2 - (672*kmrgd(1)^9*par_v4)/(kmrgd(1)^4 + par_k4^4)^3 + (384*kmrgd(1)^13*par_v4)/(kmrgd(1)^4 + par_k4^4)^4 , -(par_v2*(exp(-(kmrgd(1) - par_thetah)/par_kh)/(par_kh^2*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) - (2*exp(-(2*(kmrgd(1) - par_thetah))/par_kh))/(par_kh^2*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^3) + exp(-(kmrgd(1) - par_thetam)/par_km)/(par_km^2*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2) - (2*exp(-(2*(kmrgd(1) - par_thetam))/par_km))/(par_km^2*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^3)))/par_lambda ; (312*kmrgd(1)^5*par_v4)/(kmrgd(1)^4 + par_k4^4)^2 - (24*kmrgd(1)*par_v4)/(kmrgd(1)^4 + par_k4^4) - (672*kmrgd(1)^9*par_v4)/(kmrgd(1)^4 + par_k4^4)^3 + (384*kmrgd(1)^13*par_v4)/(kmrgd(1)^4 + par_k4^4)^4 , 0 ];
tens32=[ -(par_v2*(exp(-(kmrgd(1) - par_thetah)/par_kh)/(par_kh^2*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) - (2*exp(-(2*(kmrgd(1) - par_thetah))/par_kh))/(par_kh^2*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^3) + exp(-(kmrgd(1) - par_thetam)/par_km)/(par_km^2*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2) - (2*exp(-(2*(kmrgd(1) - par_thetam))/par_km))/(par_km^2*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^3)))/par_lambda , 0 ; 0 , 0 ];
tens33=[ -(par_v2*(exp(-(kmrgd(1) - par_thetah)/par_kh)/(par_kh^2*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^2) - (2*exp(-(2*(kmrgd(1) - par_thetah))/par_kh))/(par_kh^2*(exp(-(kmrgd(1) - par_thetah)/par_kh) + 1)^3) + exp(-(kmrgd(1) - par_thetam)/par_km)/(par_km^2*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^2) - (2*exp(-(2*(kmrgd(1) - par_thetam))/par_km))/(par_km^2*(exp(-(kmrgd(1) - par_thetam)/par_km) + 1)^3)))/par_lambda , 0 ; 0 , 0 ];
tens34=[ 0 , 0 ; 0 , 0 ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,2,1) =tens33;
tens3(:,:,2,2) =tens34;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_v1,par_v2,par_v3,par_v4,par_k3,par_k4,par_thetam,par_km,par_thetah,par_kh,par_lambda,par_j0,par_j1,par_a)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_v1,par_v2,par_v3,par_v4,par_k3,par_k4,par_thetam,par_km,par_thetah,par_kh,par_lambda,par_j0,par_j1,par_a)
