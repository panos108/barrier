clc
clear
close all
 load data3
 x=[-.9:0.005:2];%x(2,1:500);
 u=sinh((x));%u(2:501);Auto trexei lala
% x1=x;
% x=u;
% u=x1;

% for i=1:length(x)
%     if x(i)>=0
%         u(i)=sqrt(x1(i));
%     else
%         u(i)=-sqrt(x1(i));
%     end
% end
        % u(100)=0;
% x(100)=0;
for i=1:length(x)
    b(i)=0;
    db(i)=0;
    xx=x(i);
    U=u(i);
    if (xx+.1*U(1)+2)<1e-2
        b(i)=+0.5*(((xx+.1*U(1)+2-2*0.01)/0.01)^2-1)-log(0.01)+.5*(xx+.1*U(1))+log(2);
        %                 b=b+0.5*(((xx(1)+2-2*0.000001)/0.000001)^2-1)-log(0.000001)+.5*(xx(1))+log(2);
    else
        b(i)=-log(xx+.1*U(1)+2)+log(2)+0.5*(xx+.1*U(1));%+log(WW(j));%+log(WW(j)-l(j))+G(j,i)*U(i);
        %           b=b-log(xx(1)+2)+log(2)+0.5*(xx(1));%+log(WW(j));%+log(WW(j)-l(j))+G(j,i)*U(i);
        
    end
    
    xr=xx+1e-5;
    xl=xx-1e-5;
        br(i)=-log(xr+.1*U(1)+2)+log(2)+0.5*(xr+.1*U(1));%+log(WW(j));%+log(WW(j)-l(j))+G(j,i)*U(i);
        bl(i)=-log(xl+.1*U(1)+2)+log(2)+0.5*(xl+.1*U(1));%+log(WW(j));%+log(WW(j)-l(j))+G(j,i)*U(i);
        db1(i)=.5*(br(i)-bl(i))/1e-5;
    
    
    if (xx+.1*U(1)+2)<1e-2
        db(i)=+0.5*2*0.1*(xx+0.1*U+2-2*0.01)/0.01^2+0.5*0.1;
        %                 b=b+0.5*(((xx(1)+2-2*0.000001)/0.000001)^2-1)-log(0.000001)+.5*(xx(1))+log(2);
    else
        db(i)=-0.1/(xx+0.1*U+2)+0.5*0.1;%+log(WW(j));%+log(WW(j)-l(j))+G(j,i)*U(i);
        %           b=b-log(xx(1)+2)+log(2)+0.5*(xx(1));%+log(WW(j));%+log(WW(j)-l(j))+G(j,i)*U(i);
        
    end
    ud(i)=(u(i))*db(i);
    ll(i)=xx+.1*U(1)+2;
end
% plot(u,b,'+')
% hold on
% figure
% plot(u,ud,'o')
% figure
% plot(u,db1)
% % 
% % 
%  for i=2:length(x)
%      
% %     d(i-1)=(u(i)-u(i-1))*(db(i)-db(i-1));
%        if x(i)>-0.7
%     ud1(i)=u(i)*db(i)*1/[-1.64];
%        else
%     ud1(i)=u(i)*db(i)*1/[(u(i)-u(i-1))/(x(i)-x(i-1))];
%        end
% l(i)=[(u(i)-u(i-1))/(x(i)-x(i-1))];
% end
% figure
%  plot(u,ud1,'+')
% 
%  
%  k=0;
%  for i=1:500
%      k=k+1;
%      if x(i)>-0.7
%          lx(k+1)=x(i);
%          ly(k+1)=u(i);
%      end
%  end
  k = convhull(u,b);
plot(u(k),b(k),'r-',u,b,'--')


figure 
k = convhull(-x,u);
plot(-x(k),u(k),'r-',-x,u,'--')

figure 
k = convhull(u,ll);
plot(u(k),ll(k),'r-',u,ll,'--')


syms x
f=-log(asinh(x)+.1*x+2)+log(2)+0.5*(asinh(x)+.1*x);
df2=diff(f,x,2);
df=diff(f,x);
f1=tanh(x)+x;
df1=diff(f1,x);
df12=diff(f1,x,2);
l=1/(f1+2)^2*df1^2+(.5-1/(f1+2))*df12;
l1=1/(f1+2)^2*df1^2;
l2=(.5-1/(f1+2));
l3=df12;
figure
ezplot(l2);hold on;ezplot(l3)
ezplot(l);hold on;
figure
ezplot(df2,[min(u),max(u)])
figure
ezplot(l2,[min(u),max(u)]);hold on;ezplot(l3,[min(u),max(u)])





