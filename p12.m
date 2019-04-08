% p12.m - accuracy of Chebyshev spectral differentiation
% 	(compare p7.m)
% Compute derivatives for various values of N:

tic;
Nmax = 500; E = zeros(3,Nmax);
for N = 1:Nmax;
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
end
toc;

% Plot results:
titles = {'|x^3|','exp(-x^{-2})','1/(1+x^2)','x^{10}'}; clf
for iplot = 1:4
    subplot(2,2,iplot)
    semilogy(1:Nmax,E(iplot,:),'.','markersize',12)
    line(1:Nmax,E(iplot,:),'linewidth',.8)
    axis([0 Nmax 1e-16 1e3]), grid on
    set(gca,'xtick',0:10:Nmax,'ytick',(10).^(-15:5:0))
    xlabel N, ylabel error, title(titles(iplot))
end
