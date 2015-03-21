function val = lagrangian(alpha, beta, gamma, x, y, z, X, Y, S, T, rho, n, p, d)

lag = sum(gamma);
for i = 1:n,
    temp = zeros(2,2);
    for j = 1:p,
        temp = temp + alpha(j,1)*[0 0; 0 -z(j,i)];
    end

    for j = 1:d,
        temp = temp + beta(j,1)*[0 x(j,i); x(j,i) 0];
    end

    temp = temp + gamma(i,1)*[-1 0; 0 0];

    temp = temp + S(:,:,i);
    temp = temp - [0 y(1,i); y(1,i) 0];

    lag = lag - trace(X(:,:,i)'*temp) + (rho/2)*norm(temp,'fro')^2;
end

temp = zeros(p+1,1);
temp(1:p,1) = -alpha + T(1:p,1);
temp(p+1,1) = sum(alpha) + T(p+1,1) - 1;
    
val = lag - sum(Y.*temp) + (rho/2)*norm(temp,'fro')^2';

end