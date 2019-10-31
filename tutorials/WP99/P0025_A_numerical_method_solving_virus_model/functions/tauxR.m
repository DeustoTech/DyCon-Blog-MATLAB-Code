    function w=tauxR(N,h,x,delta,kk,c0,mu,sigma,alpha1)
sumr=0;
for rr=1:N
    sumr=sumr+h*feval(@p,x(rr))*exp(-x(rr)*delta);
end
w=(sumr*kk*sigma)/(c0*(mu+alpha1*sigma));
    end