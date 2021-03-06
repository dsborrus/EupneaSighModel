p=[0;0;0;3];ap1=[1];
[x0,v0]=init_EP_EP(@cstr,[-0.9],p,ap1);
opt=contset;
opt=contset(opt,'VarTolerance',1e-3);
opt=contset(opt,'FunTolerance',1e-3);
opt=contset(opt,'MaxNumPoints',50);
opt=contset(opt,'Singularities',1);
[x,v,s,h,f]=cont(@equilibrium,x0,[],opt);
cpl(x,v,s,[2,1]);
x1=x(1,s(3).index);
p(ap1)=x(end,s(3).index);
[x0,v0]=init_LP_LP(@cstr,x1,p,[1 2],[1 2 3 4]);
opt=contset(opt,'MaxNumPoints',300);
[x2,v2,s2,h2,f2]=cont(@limitpoint,x0,v0,opt);
hold on;
cpl(x2,v2,s2,[2,1]);
x1=x2(1,s2(2).index);
p([1 2])=x2(end-1:end,s2(2).index);
[x0,v0]=init_BP_BP(@cstr,x1,p,[1 2 3],1);
opt=contset(opt,'Backward',1);
[x3,v3,s3,h3,f3]=cont(@branchpoint,x0,[],opt);
hold on;
cpl(x3,v3,s3,[2,1]);
x1=x2(1,s2(5).index);
p([1 2])=x2(end-1:end,s2(5).index);
[x0,v0]=init_BP_BP(@cstr,x1,p,[1 2 3],4);
opt=contset(opt,'Backward',1);
[x3,v3,s3,h3,f3]=cont(@branchpoint,x0,[],opt);
hold on;
cpl(x3,v3,s3,[2,1]);