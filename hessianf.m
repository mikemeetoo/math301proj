function [h]=hessianf(alpha,beta,c);
    [p,mm]=size(alpha);
    [d,mm]=size(beta);
    [n,mm]=size(c);
    X=c(:,1:d);
    Y=c(:,d+1);
    Z=c(:,d+2:d+1+p);
    w1=ones(n,1)./(Z*alpha);
    w2=X*beta-Y;
    w = w2.*w1;
    WW=sparse(1:n,1:n,w1.*(w.^2));
    h11=Z'*WW*Z;
    WW=sparse(1:n,1:n,w1.*w);
    h12=-Z'*WW*X;
    WW=sparse(1:n,1:n,w1);
    h22= X'*WW*X;
    h=2*[h11 h12;h12' h22];
return
