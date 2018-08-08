% Els Fieremans, Hong-Hsi Lee, Physical and numerical phantoms for the
% validation of brain microstructural MRI: A cookbook, NeuroImage 2018
%
% Supplementary Information: Recipes of MC simulations in Figure 4
%
% Example 1: Free diffusion in 2d (Figure 4, point 1)

clear

% Setup simulation parameters
D0 = 2;             % Intrinsic diffusivity
dt = 5e-2;          % Time for each step
Ns = 2e3;           % # Steps
Np = 5e4;           % # Particles
t = (1:Ns)*dt;      % Diffusion time
dx = sqrt(4*D0*dt); % Step size

% Diffusion phase
phi = rand(Np,Ns)*2*pi;

% Diffusion displacement and cumulant in x-direction
x = cos(phi);
x = cumsum(x,2);
x2 = sum(x.^2,1)*dx^2/Np;
x4 = sum(x.^4,1)*dx^4/Np;

% Diffusivity and Kurtosis in x-direction
Dx = x2/2./t;
Kx = x4./x2.^2-3;

%% Plot figures

% Plot time-dependent diffusivity in x-direction
figure; subplot(131)
h = plot(t,Dx/D0); set(h,'linewidth',3)
xlim([0 max(t)]); ylim([0 1.25])
box on; grid on
set(gca,'fontsize',20,'xtick',0:20:100,'ytick',0:0.5:3)
pbaspect([1 1 1])
xlabel('$t$ (ms)','interpreter','latex','fontsize',30)
ylabel('$D(t)/D_0$','interpreter','latex','fontsize',30)

% Plot time-dependent kurtosis in x-direction
subplot(132);
h = plot(t,Kx); set(h,'linewidth',3)
xlim([0 max(t)]); ylim([-2 0.25])
box on; grid on
set(gca,'fontsize',20,'xtick',0:20:100,'ytick',-3:0.5:3)
pbaspect([1 1 1])
xlabel('$t$ (ms)','interpreter','latex','fontsize',30)
ylabel('$K(t)$','interpreter','latex','fontsize',30)

% Plot time-dependent kutosis wrt the step number to demonstrate K~1/Ns !
subplot(133);
h = plot(1./(1:Ns),Kx,'.'); set(h,'markersize',20)
hold on;
hr = refline(-1.5,0); set(hr,'linewidth',1,'color','r')
legend([h,hr],{'Simulation','Theory'},'interpreter','latex','fontsize',20)
xlim([0 1]); ylim([-2 0.25])
box on; grid on
set(gca,'fontsize',20,'xtick',0:0.2:1,'ytick',-3:0.5:3)
pbaspect([1 1 1])
xlabel('$1/N_{\rm{step}}$','interpreter','latex','fontsize',30)
ylabel('$K(t)$','interpreter','latex','fontsize',30)

set(gcf,'unit','inch','position',[0 0 15 5])

