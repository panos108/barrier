clear all
clc
x=[-1:0.01:1];
global h x1 mi1
h=0.5;
for i=1:length(x)
u(i)=quadprog(h,-x(i),[],[],[],[],-1,1);
end
plot(x,u)
mi(1)=1;
for j=1:100
    mi1=mi(j);
for i=1:length(x)
    x1=x(i);
[U,J,F]=fminsearch(@funsimple,zeros(1,1));
u1(i)=U;
uu(i,j)=u1(i);
end
s1(j)=(u1(end)-u1(1))/(x(end)-x(1));
h1(j)=1/s1(j);
if j<100
mi(j+1)=mi(j)-0.01;
end
end
hold on
plot(x,uu')
[p,S]=polyfit(mi,h1,1);
yfit =  p(1) * mi + p(2);
yresid = h1 - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(h1)-1) * var(h1);
rsq = 1 - SSresid/SStotal