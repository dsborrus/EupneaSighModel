jin0=0.02;
jin1=0;
v4=0.3;
n4=4;
k4=0.3;
v1=20;

figure(5)
cc = 0:0.001:0.4;
ctct = 0:0.001:3;
css0=fzero(@(x) jin0-v4*x^n4/(k4^n4+x^n4),0.5);
plot(css0*ones(size(ctct)),ctct,'k'); hold on;
css1=fzero(@(x) jin0+jin1-v4*x^n4/(k4^n4+x^n4),0.5);
plot(css1*ones(size(ctct)),ctct,'y-.'); hold on;
out.cnull_ct=ctct;
out.cnull_c=css0*ones(size(ctct));

for i=1:length(cc)
    ctnull(i)=fzero(@(x) (v2+v1*finf(cc(i),thetam,km,thetah,kh))*((x-cc(i))/lambda-cc(i))-v3*cc(i)^n3/(k3^n3+cc(i)^n3),0.3);
end
plot(cc,ctnull,'r'); hold on;
plot(c(:),ct(:),'b-o'); hold on; xlabel('c'); ylabel('ct'); title(titlestr)
plot(cc,max(ctct)*cc.^n4./(k4^n4+cc.^n4),'c');
print([ run '-5.png'],'-dpng')
out.ctnull_c=cc;
out.ctnull_ct=ctnull;