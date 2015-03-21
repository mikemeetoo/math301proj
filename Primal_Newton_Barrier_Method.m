clear
clc 

randn('seed', 300);

%Test with randomly generated data
%{

p = 100;
d = 100;
n = 400;

c = zeros(n,p+d+1);

c(:,1:d) = randn(n,d);
c(:,d+1) = randn(n,1);
c(:,d+2:end) = abs(randn(n,p));

xh = c(:,1:d)';
yh = c(:,d+1)';
zh = c(:,d+2:end)';

tic
cvx_begin 
        variable alp(p) nonnegative
        variables bet(d) gam(n,1);
        
        sum(alp) == 1;
        
        for i = 1:n,
            [gam(i,1) yh(1,i)-xh(:,i)'*bet;
             yh(1,i)-xh(:,i)'*bet zh(:,i)'*alp] == semidefinite(2);
        end
        minimize ( sum(gam) )
            
cvx_end
toc 

cvx_optval

x0 = ones(p,1)./p;
z0 = zeros(d,1);
y0 = -4000;
%}

load regret_data_owen.mat
p = 51;
d = 50;
n = 200;

%Initialize
A = ones(1,p);
b = 1;

opt_val_store = [];

tol=1.e-6; 
gamma=.6;    
alpha=0.3;    

%Primal Newton Barrier Method
x=x0;
z=z0;
y=y0;
g=gradienf(x,z,c);

s=g(1:p)-A'*y;
mu = x'*s/p;
iter =0; 

tic
while mu >= tol,
    iter = iter + 1;
    
    mu = gamma*mu;
    
    xr = ones(p,1)./x;
    lhs = [mu*xr-s;-g(p+1:p+d);-(A*x-b)]; %Ax-b can be change to 0

    Hess = hessianf(x,z,c);
   
    Hess(1:p,1:p)= Hess(1:p,1:p)+mu*sparse(1:p,1:p,xr.^2);
 
    sol = [Hess [-A';sparse(d,1)];A sparse(1,d+1)]\lhs;
  
    dx= sol(1:p);
    dz= sol(p+1:p+d);
    dy= sol(p+d+1:p+d+1);
  
    step_size = 1;
    ss = 0;
    
    while ss <= 0,
        xx = x + step_size*dx;
        yy = y + step_size*dy;
        zz = z + step_size*dz;
        g=gradienf(xx,zz,c);
        s = g(1:p) - A'*yy;
        step_size = alpha*step_size;
        ss = min([xx;s]);
    end
    
    x=xx;
    y=yy;
    z=zz;
    mu = x'*s/p;
    
    %Probe Objective Value
    %{
    alphah = x;
    betah = z;

    xh = c(:,1:d)';
    yh = c(:,d+1)';
    zh = c(:,d+2:end)';

    sm = 0;
    for i = 1:n,
     sm = sm + (yh(1,i)-xh(:,i)'*betah)^2/(zh(:,i)'*alphah);
    end
    opt_val = sm;
    opt_val_store = [opt_val_store norm(opt_val-cvx_optval)];
    %}
end
toc

alphah = x;
betah = z;

xh = c(:,1:d)';
yh = c(:,d+1)';
zh = c(:,d+2:end)';

iter
sm = 0;
for i = 1:n,
 sm = sm + (yh(1,i)-xh(:,i)'*betah)^2/(zh(:,i)'*alphah);
end
opt_val = sm


%{
figure();
hold on;
plot(1:iter, opt_val_store);
hXLabel = xlabel('Iter');
hLegend = ylabel('||f(x^*) - f(x_k)||');
hTitle = title('Pimal Newton Barrier Method');
%}



