function z=int(z,N,h,x)
% Calcuiate the integral
sum=0;
for l=1:N+1
    sum=sum+h*feval(@p,x(l))*z(l);
end
    z=sum;
end