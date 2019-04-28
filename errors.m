Nmax = 50;

e1 = zeros(Nmax, 1);
e2 = zeros(Nmax, 1);
e3 = zeros(Nmax, 1);

% Chebychev differentiation via matrices
for N = 1:Nmax;
    [D,x] = cheb(N);
    v = exp(-x.^(-2)); vprime = 2.*v./x.^3; % C-infinity
    e1(N) = norm(D*v-vprime,inf);
end

% Chebychev differentiation via FFT
for N = 1:Nmax;
    x = cos(pi*(0:N)'/N);
    v = exp(-x.^(-2)); vprime = 2.*v./x.^3; % C-infinity
    e2(N) = norm(chebfft(v, x)-vprime,inf);
end

for N = 2:Nmax;
    h = 2*pi/N; x = -pi + (1:N)'*h;
    % Construct sparse 4th-order differentiation matrix:
    e = ones(N,1);
    D = sparse(1:N,[2:N 1],2*e/3,N,N) - sparse(1:N,[3:N 1 2],e/12,N,N);
    D = (D-D')/h;

    v = exp(-x.^(-2)); vprime = 2.*v./x.^3; % C-infinity
    e3(N) = norm(D*v-vprime,inf);
end

clf;
plot(1:Nmax, e1)
hold on
plot(1:Nmax, e2)
plot(1:Nmax, e3, '*b')
hold off
title('Rate of convergence')
xlabel('N')
ylabel('error')
legend('matrix','FFT','finite difference')
