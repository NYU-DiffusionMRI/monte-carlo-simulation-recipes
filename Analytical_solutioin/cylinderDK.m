function [RD,RK] = cylinderDK(r,D0,t,N)
%CYLINDERDK    radial diffusivity and radial kurtosis inside a cylinder
%   [RD,RK] = CYLINDERDK(r,D0,t,N) returns radial diffusivity RD and radial
%   kurtosis RK inside a cylinder in radius r at diffusion time t. The
%   intrinsic diffusivity inside the cylinder is D0, and the number of
%   summation terms is N.
%
%   Author: Hong-Hsi Lee, September, 2018 (orcid.org/0000-0002-3663-6559)

dJ0 = @(x) -besselj(1,x);
beta0k = @(k) fzero(dJ0, [(k-1)*pi, k*pi]); 
beta0 = zeros(1,N); for k=1:N, beta0(k) = beta0k(k+1); end
% beta0 = [3.8317    7.0156   10.1735   13.3237   16.4706   19.6159   22.7601   25.9037   29.0468   32.1897];
beta0 = [NaN beta0(1:N-1)];

dJ1 = @(x) besselj(0,x) - besselj(1,x)./(x + eps);
beta1k = @(k) fzero(dJ1, [(k-1)*pi, k*pi]); 
beta1 = zeros(1,N); for k=1:N, beta1(k) = beta1k(k); end
% beta1 = [1.8412    5.3314    8.5363   11.7060   14.8636   18.0155   21.1644   24.3113   27.4571   30.6019];

dJ2 = @(x) besselj(1,x) - besselj(2,x)*2./(x + eps);
beta2k = @(k) fzero(dJ2,[floor((k-1)*pi), k*pi]);
beta2 = zeros(1,N); for k=1:N, beta2(k) = beta2k(k+1); end
% beta2 = [3.0542    6.7061    9.9695   13.1704   16.3475   19.5129   22.6716   25.8260   28.9777   32.1273];
beta2 = [NaN beta2(1:N-1)];

RD = r^2/4./t;
tc = r^2/D0;
for i = 1:numel(beta1)
    RD = RD - 2*r^2./t.*exp(-beta1(i)^2*t/tc)/beta1(i)^2/(beta1(i)^2-1);
end

numer = 5/64/3 + exp(-beta1(1)^2*t/tc)/beta1(1)^2/(beta1(1)^2-1)*(4/beta1(1)^2-1.5);
denom = RD./(r^2./t);
for i = 2:numel(beta0)
    numer = numer + exp(-beta0(i)^2*t/tc)/beta0(i)^4 +...
        exp(-beta1(i)^2*t/tc)/beta1(i)^2/(beta1(i)^2-1)*(4/beta1(i)^2-3/2) +...
        exp(-beta2(i)^2*t/tc)/beta2(i)^2/(beta2(i)^2-4)/2;
end
numer = 6*numer;
denom = (denom).^2;
RK = numer./denom-3;

end