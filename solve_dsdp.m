load regret_data_owen.mat
p = 51;
d = 50;
n = 200;

X = c(:,1:d);
Y = c(:,d+1);
Z = c(:,d+2:end);

for i=1:n
    
    %Input 2x2 matrix coefficients of alpha
    AACi=sparse(0,0);
    AACi = [ AACi [sparse(2,p); -sparse(Z(i,:))] ];
    
    %Input 2x2 matrix coefficients of beta
    AACi = [ AACi [sparse(1,d); sparse(X(i,:)); sparse(1,d)] ];
    
    %Input 2x2 matrix coefficients of gamma(i)
    temp = sparse(3,n);
    temp(1,i) = -1;
    AACi = [ AACi temp ];
    
    %Input 2x2 matrix of the ith right-hand-side of the dual
    AACi=[ AACi sparse([ 0 Y(i) 0 ]') ];
    
    AC{i,1} = 'SDP';
    AC{i,2} = [2];
    AC{i,3} = AACi;
end

%Set up nonnegative cone for alpha
AC{n+1,1}='LP';
AC{n+1,2}=p+1;
AC{n+1,3}=[[sparse(ones(1,p));-speye(p)] sparse(p+1,d) sparse(p+1,n) sparse([1 sparse(1,p)]')];

%Set up the dual objective vector
b=[zeros(1,p) zeros(1,d) -ones(1,n)]';

[STAT,y,XX]=dsdp(b,AC);
alpha=y(1:p);
beta =y(p+1:p+d);
