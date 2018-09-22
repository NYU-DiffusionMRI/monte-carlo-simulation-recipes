function [D,K] = parallelplaneDK(d,D0,t,N)
%PARALLELPLANEDK    diffusivity and kurtosis between parallel planes
%   [D,K] = PARALLELPLANEDK(d,D0,t,N) returns diffusivity D and kurtosis K
%   between two parallel planes with a spacing d at diffusion time t. The
%   intrinsic diffusivity between the planes is D0, and the number of
%   summation terms is N.
%
%   Author: Hong-Hsi Lee, September, 2018 (orcid.org/0000-0002-3663-6559)

a = d/2;

dXi = @(x) x.*tan(x);
xin = @(n) fzero(dXi, [(n-1/2)*pi+1e-10, (n+1/2)*pi-1e-10]);
xi = zeros(1,N); for n = 1:N, xi(n) = xin(n); end
% xi = pi*(1:10);
xi = [NaN, xi(1:N-1)];

dZi = @(x) x.*cot(x);
zetam = @(m) fzero(dZi, [(m-1)*pi+1e-10, m*pi-1e-10]);
% zi = pi*(0.5:9.5);
zeta = zeros(1,N); for m = 1:N, zeta(m) = zetam(m); end

D = a^2/3./t;
tc = a^2/D0;
for i = 1:numel(zeta)
    D = D - 2*a^2./t.*exp(-zeta(i)^2*t/tc)/zeta(i)^4;
end

numer = 2/45 + exp(-zeta(1)^2*t/tc)*(2/zeta(1)^4)*(-1+2/zeta(1)^2);
denom = D./(a^2./t);
for i = 2:numel(xi)
    numer = numer + exp(-xi(i)^2*t/tc)*2/xi(i)^4 +...
        exp(-zeta(i)^2*t/tc)*(2/zeta(i)^4)*(-1+2/zeta(i)^2);
end
numer = 6*numer;
denom = (denom).^2;
K = numer./denom-3;

end