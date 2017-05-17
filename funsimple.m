function J=funsimple(U)
global h x1 mi1

% J=0;
% K=WW-G*U;
% l=G*U;
    b=0;
for i=1:1
if U(i)+1>0.0000001
        b=b-log((U(i)+2));%-0.5*u(j);       
elseif    U(i)+1<0.0000001
        b=b+0.5*(((U(i)+1-2*0.0000001)/0.0000001)^2-1)-log(0.0000001);
end
if -U(i)+1>0.0000001
        b=b-log((-U(i)+1));%-0.5*u(j);       
elseif    -U(i)+1<0.000001
        b=b+0.5*(((-U(i)+1-2*0.0000001)/0.0000001)^2-1)-log(0.0000001);
end
% -log((1-U(i)))
%     b1=0;
%     for j=1:1
%         b=b-log((1-U(i)))-log((U(i)+1));%-0.5*u(j);       
%     end

%     b1=b1+log(3)-log(3-x(1,i))+1/3*x(1,i);
%    J=J+.5*x(:,i+1)'*x(:,i+1)+.5*u(i)'*0.1*u(i)+1e-10*b;
end


J=.5*U'*h*U-x1'*U+mi1*b;