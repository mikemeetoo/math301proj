function [g]=gradienf(alpha,beta,c);
    [p,mm]=size(alpha);
    [d,mm]=size(beta);
    [n,mm]=size(c);
    X=c(:,1:d);
    Y=c(:,d+1);
    Z=c(:,d+2:d+1+p);
    w1=ones(n,1)./(Z*alpha);
    w2=X*beta-Y;
    w=w2.*w1;
    g1=-Z'*(w.^2);
    g2= 2*(X'*w);
    g=[g1;g2];
return
