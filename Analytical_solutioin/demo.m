% Els Fieremans, Hong-Hsi Lee, Physical and numerical phantoms for the
% validation of brain microstructural MRI: A cookbook, NeuroImage 2018
%
% Appendix E. Analytical solutions for the diffusion time-dependence in
% specific geometries
%
% 1. Between two parallel planes
% 2. Inside a cylinder
% 3. Inside a sphere

%% Between two parallel planes

d = 2;      % spacing between planes, µm
D0 = 2;     % intrinsic diffusivity, µm^2/ms
N = 100;    % # summation term

% short-time limit
t1 = 1e-4:1e-3:2;
[D1,K1] = parallelplaneDK(d,D0,t1,N);
SV = 2/d;   % surface-to-volume ratio
D_mitra = D0*(1-SV/1*4/3/sqrt(pi)*sqrt(D0*t1));
K_mitra = SV*8/5/sqrt(pi)*sqrt(D0*t1);

% long-time limit
t2 = 1:0.05:100;
[D2,K2] = parallelplaneDK(d,D0,t2,N);
D_longtime = d^2./t2/12;
K_longtime = -3/5*ones(size(t2));

% plot figure
figure;
subplot(221); hold on
h = plot(sqrt(t1),D1/D0,'.');
h_mitra = plot(sqrt(t1),D_mitra/D0,'--');
set(h_mitra,'linewidth',1);
box on; grid on;
xlim([0,sqrt(max(t1))]); ylim([0,1]);
legend([h,h_mitra],{'Simulation','Short-time limit'},'interpreter','latex','fontsize',12)
xlabel('$\sqrt{t}$ (ms$^{1/2}$)','interpreter','latex','fontsize',12)
ylabel('$D(t)/D_0$','interpreter','latex','fontsize',12)
title('Short-time limit','interpreter','latex','fontsize',12)

subplot(222); hold on
h = plot(sqrt(t1),K1,'.');
h_mitra = plot(sqrt(t1),K_mitra,'--');
set(h_mitra,'linewidth',1);
box on; grid on
xlim([0,sqrt(max(t1))]); ylim([-1,0.5]);
legend([h,h_mitra],{'Simulation','Short-time limit'},'interpreter','latex','fontsize',12)
xlabel('$\sqrt{t}$ (ms$^{1/2}$)','interpreter','latex','fontsize',12)
ylabel('$K(t)$','interpreter','latex','fontsize',12)
title('Short-time limit','interpreter','latex','fontsize',12)

subplot(223); hold on
h = plot(1./t2,D2/D0,'.');
h_longtime = plot(1./t2,D_longtime/D0,'--');
set(h_longtime,'linewidth',1);
box on; grid on;
xlim([0,1]); ylim([0,0.2]);
legend([h,h_longtime],{'Simulation','Long-time limit'},'interpreter','latex','fontsize',12)
xlabel('$1/t$ (ms$^{-1}$)','interpreter','latex','fontsize',12)
ylabel('$D(t)/D_0$','interpreter','latex','fontsize',12)
title('Long-time limit','interpreter','latex','fontsize',12)

subplot(224); hold on
h = plot(t2,K2,'.');
h_longtime = plot(t2,K_longtime,'--');
set(h_longtime,'linewidth',1);
box on; grid on
xlim([0,100]); ylim([-1,1]);
legend([h,h_longtime],{'Simulation','Long-time limit'},'interpreter','latex','fontsize',12)
xlabel('$t$ (ms)','interpreter','latex','fontsize',12)
ylabel('$K(t)$','interpreter','latex','fontsize',12)
title('Long-time limit','interpreter','latex','fontsize',12)

%% Inside a cylinder

r = 1;      % cylinder radius, µm
D0 = 2;     % intrinsic diffusivity, µm^2/ms
N = 100;    % # summation term

% short-time limit
t1 = 1e-4:1e-3:2;
[D1,K1] = cylinderDK(r,D0,t1,N);
SV = 2/r;   % surface-to-volume ratio
D_mitra = D0*(1-SV/2*4/3/sqrt(pi)*sqrt(D0*t1) - SV/4/2/r*D0*t1);
K_mitra = SV*(3/8)*8/5/sqrt(pi)*sqrt(D0*t1);

% long-time limit
t2 = 1:0.05:100;
[D2,K2] = cylinderDK(r,D0,t2,N);
D_longtime = r^2./t2/4;
K_longtime = -1/2*ones(size(t2));

% plot figure
figure;
subplot(221); hold on
h = plot(sqrt(t1),D1/D0,'.');
h_mitra = plot(sqrt(t1),D_mitra/D0,'--');
set(h_mitra,'linewidth',1);
box on; grid on;
xlim([0,sqrt(max(t1))]); ylim([0,1]);
legend([h,h_mitra],{'Simulation','Short-time limit'},'interpreter','latex','fontsize',12)
xlabel('$\sqrt{t}$ (ms$^{1/2}$)','interpreter','latex','fontsize',12)
ylabel('$D(t)/D_0$','interpreter','latex','fontsize',12)
title('Short-time limit','interpreter','latex','fontsize',12)

subplot(222); hold on
h = plot(sqrt(t1),K1,'.');
h_mitra = plot(sqrt(t1),K_mitra,'--');
set(h_mitra,'linewidth',1);
box on; grid on
xlim([0,sqrt(max(t1))]); ylim([-1,0.5]);
legend([h,h_mitra],{'Simulation','Short-time limit'},'interpreter','latex','fontsize',12)
xlabel('$\sqrt{t}$ (ms$^{1/2}$)','interpreter','latex','fontsize',12)
ylabel('$K(t)$','interpreter','latex','fontsize',12)
title('Short-time limit','interpreter','latex','fontsize',12)

subplot(223); hold on
h = plot(1./t2,D2/D0,'.');
h_longtime = plot(1./t2,D_longtime/D0,'--');
set(h_longtime,'linewidth',1);
box on; grid on;
xlim([0,1]); ylim([0,0.2]);
legend([h,h_longtime],{'Simulation','Long-time limit'},'interpreter','latex','fontsize',12)
xlabel('$1/t$ (ms$^{-1}$)','interpreter','latex','fontsize',12)
ylabel('$D(t)/D_0$','interpreter','latex','fontsize',12)
title('Long-time limit','interpreter','latex','fontsize',12)

subplot(224); hold on
h = plot(t2,K2,'.');
h_longtime = plot(t2,K_longtime,'--');
set(h_longtime,'linewidth',1);
box on; grid on
xlim([0,100]); ylim([-1,1]);
legend([h,h_longtime],{'Simulation','Long-time limit'},'interpreter','latex','fontsize',12)
xlabel('$t$ (ms)','interpreter','latex','fontsize',12)
ylabel('$K(t)$','interpreter','latex','fontsize',12)
title('Long-time limit','interpreter','latex','fontsize',12)

%% Inside a sphere

r = 1;      % sphere radius, µm
D0 = 2;     % intrinsic diffusivity, µm^2/ms
N = 100;    % # summation term

% short-time limit
t1 = 1e-4:1e-3:2;
[D1,K1] = sphereDK(r,D0,t1,N);
SV = 3/r;   % surface-to-volume ratio
D_mitra = D0*(1-SV/3*4/3/sqrt(pi)*sqrt(D0*t1) - SV/4/3*(2/r)*D0*t1);
K_mitra = SV*(1/5)*8/5/sqrt(pi)*sqrt(D0*t1);

% long-time limit
t2 = 1:0.05:100;
[D2,K2] = sphereDK(r,D0,t2,N);
D_longtime = r^2./t2/5;
K_longtime = -3/7*ones(size(t2));

% plot figure
figure;
subplot(221); hold on
h = plot(sqrt(t1),D1/D0,'.');
h_mitra = plot(sqrt(t1),D_mitra/D0,'--');
set(h_mitra,'linewidth',1);
box on; grid on;
xlim([0,sqrt(max(t1))]); ylim([0,1]);
legend([h,h_mitra],{'Simulation','Short-time limit'},'interpreter','latex','fontsize',12)
xlabel('$\sqrt{t}$ (ms$^{1/2}$)','interpreter','latex','fontsize',12)
ylabel('$D(t)/D_0$','interpreter','latex','fontsize',12)
title('Short-time limit','interpreter','latex','fontsize',12)

subplot(222); hold on
h = plot(sqrt(t1),K1,'.');
h_mitra = plot(sqrt(t1),K_mitra,'--');
set(h_mitra,'linewidth',1);
box on; grid on
xlim([0,sqrt(max(t1))]); ylim([-1,0.5]);
legend([h,h_mitra],{'Simulation','Short-time limit'},'interpreter','latex','fontsize',12)
xlabel('$\sqrt{t}$ (ms$^{1/2}$)','interpreter','latex','fontsize',12)
ylabel('$K(t)$','interpreter','latex','fontsize',12)
title('Short-time limit','interpreter','latex','fontsize',12)

subplot(223); hold on
h = plot(1./t2,D2/D0,'.');
h_longtime = plot(1./t2,D_longtime/D0,'--');
set(h_longtime,'linewidth',1);
box on; grid on;
xlim([0,1]); ylim([0,0.2]);
legend([h,h_longtime],{'Simulation','Long-time limit'},'interpreter','latex','fontsize',12)
xlabel('$1/t$ (ms$^{-1}$)','interpreter','latex','fontsize',12)
ylabel('$D(t)/D_0$','interpreter','latex','fontsize',12)
title('Long-time limit','interpreter','latex','fontsize',12)

subplot(224); hold on
h = plot(t2,K2,'.');
h_longtime = plot(t2,K_longtime,'--');
set(h_longtime,'linewidth',1);
box on; grid on
xlim([0,100]); ylim([-1,1]);
legend([h,h_longtime],{'Simulation','Long-time limit'},'interpreter','latex','fontsize',12)
xlabel('$t$ (ms)','interpreter','latex','fontsize',12)
ylabel('$K(t)$','interpreter','latex','fontsize',12)
title('Long-time limit','interpreter','latex','fontsize',12)
