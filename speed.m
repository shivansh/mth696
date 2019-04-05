Nmax = 200;

t1 = zeros(Nmax, 1);
t2 = zeros(Nmax, 1);

for N = 1:Nmax;
    % Chebychev differentiation via matrices
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
    
    tStart = tic;
    
    % Chebychev differentiation via FFT
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
    
    tElapsed = toc(tStart);
    t2(N) = toc(tStart);  % time elapsed
end

plot(1:Nmax, t1)
hold on
plot(1:Nmax, t2)
hold off
legend('matrix','FFT')