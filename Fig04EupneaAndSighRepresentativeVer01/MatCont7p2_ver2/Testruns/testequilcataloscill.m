p=[2.5;2.204678;10;0.0675;1;0.1;0.4];
ap1=[2];
[x0,v0]=init_EP_EP(@cataloscill,[0.001137;0.891483;0.062345],p,ap1);
opt=contset;
opt=contset(opt,'MaxStepSize',0.025);
opt=contset(opt,'MaxNumPoints',78);
opt=contset(opt,'Singularities',1);
[x,v,s,h,f]=cont(@equilibrium,x0,[],opt);
cpl(x,v,s,[4,1]);
