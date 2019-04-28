Nmax = 500;

t1 = zeros(Nmax, 1);
t2 = zeros(Nmax, 1);
t3 = zeros(Nmax, 1);

% Chebychev differentiation via matrices
for N = 1:Nmax;
    tStart = tic;

    [D,x] = cheb(N);
    v = abs(x).^3; vprime = 3*x.*abs(x);
    % 3rd deriv in BV
    E(1,N) = norm(D*v-vprime,inf);
    v = exp(-x.^(-2)); vprime = 2.*v./x.^3; % C-infinity
    E(2,N) = norm(D*v-vprime,inf);
    v = 1./(1+x.^2); vprime = -2*x.*v.^2;
    % analytic in [-1,1]
    E(3,N) = norm(D*v-vprime,inf);
    v = x.^10; vprime = 10*x.^9;
    % polynomial
    E(4,N) = norm(D*v-vprime,inf);

    t1(N) = toc(tStart);  % time elapsed
end

% Chebychev differentiation via FFT
for N = 1:Nmax;
    tStart = tic;

    x = cos(pi*(0:N)'/N);
    v = abs(x).^3; vprime = 3*x.*abs(x);
    % 3rd deriv in BV
    E(1,N) = norm(chebfft(v, x)-vprime,inf);
    v = exp(-x.^(-2)); vprime = 2.*v./x.^3; % C-infinity
    E(2,N) = norm(chebfft(v, x)-vprime,inf);
    v = 1./(1+x.^2); vprime = -2*x.*v.^2;
    % analytic in [-1,1]
    E(3,N) = norm(chebfft(v, x)-vprime,inf);
    v = x.^10; vprime = 10*x.^9;
    % polynomial
    E(4,N) = norm(chebfft(v, x)-vprime,inf);

    t2(N) = toc(tStart);  % time elapsed
end

for N = 2:Nmax;
    tStart = tic;

    h = 2*pi/N; x = -pi + (1:N)'*h;
    % Construct sparse 4th-order differentiation matrix:
    e = ones(N,1);
    D = sparse(1:N,[2:N 1],2*e/3,N,N) - sparse(1:N,[3:N 1 2],e/12,N,N);
    D = (D-D')/h;

    v = abs(x).^3; vprime = 3*x.*abs(x);
    err1 = norm(D*v-vprime,inf);

    v = exp(-x.^(-2)); vprime = 2.*v./x.^3; % C-infinity
    err2 = norm(D*v-vprime,inf);

    v = 1./(1+x.^2); vprime = -2*x.*v.^2;
    err3 = norm(D*v-vprime,inf);

    v = x.^10; vprime = 10*x.^9;
    err4 = norm(D*v-vprime,inf);

    t3(N) = toc(tStart);  % time elapsed
end

plot(1:Nmax, t1)
hold on
plot(1:Nmax, t2)
plot(1:Nmax, t3)
hold off
title('Runtime comparison')
xlabel('N')
ylabel('time')
legend('matrix','FFT','finite difference')
