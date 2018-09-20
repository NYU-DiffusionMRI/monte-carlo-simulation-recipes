% Els Fieremans, Hong-Hsi Lee, Physical and numerical phantoms for the
% validation of brain microstructural MRI: A cookbook, NeuroImage 2018
%
% Supplementary Information: Recipes of MC simulations in Figure 4
%
% Example 3: Check against known analytical formulas for impermeable
% (non-)absorbing membranes

% ********** Setup the directory on your computer **********
root = 'your_directory_to_this_demo/example3_analytical_formula';

% Create files for simualtion
% 1. packing directory
root_packing = fullfile(root,'packing');
fileID = fopen(fullfile(root,'cpp','root.txt'),'w');
fprintf(fileID,root_packing);
fclose(fileID);

% 2. Simulation parameters
dt = 7.5e-5;  % time of each step, ms
TN = 2e4;   % # steps
Din = 2;    % IAS diffusivity, µm^2/ms
Dout = 2;   % EAS diffusivity, µm^2/ms
kappa = 0;  % permeability, µm/ms
Dmr = 0;    % myelin diffusivity, radial, µm^2/ms
Dmc = 0;    % myelin diffusivity, circumferential, µm^2/ms
Dmz = 0;    % myelin diffusivity, along z-axis, µm^2/ms
fileID = fopen(fullfile(root,'cpp','simParamInput.txt'),'w');
fprintf(fileID,sprintf('%g\n%u\n%g\n%g\n%g\n%g\n%g\n%g\n',dt,TN,Din,Dout,...
    kappa,Dmr,Dmc,Dmz));
fclose(fileID);

%% Have a look for the microstructure
% The lookup table A saves two axon labels in one integer. If the first
% and the second axon labels are ax1 and ax2, ax1 = mod(A,Nmax), and ax2 =
% floor(A/Nmax).
% Other parameters:
%   fov: field of view of the entire gemoetry and lookup table A in µm
%   Nax: # axons
%   rCir: axon's outer radius
%   gratio: g-ratio, the ratio of inner to outer radius
%   [xCir,yCir]: axon's center position

A = load(fullfile(root_packing,'phantom_APix.txt'));
Nmax = load(fullfile(root_packing,'phantom_Nmax.txt'));
fov = load(fullfile(root_packing,'phantom_res.txt'));

Nax = load(fullfile(root_packing,'phantom_NAx.txt'));
gratio = load(fullfile(root_packing,'phantom_gratio.txt')); 
rCir = load(fullfile(root_packing,'phantom_rCir.txt'));
xCir = load(fullfile(root_packing,'phantom_xCir.txt'));
yCir = load(fullfile(root_packing,'phantom_yCir.txt'));

%% Plot the microstructure and the lookup table
% Plot the microstructure
figure; set(gcf,'unit','inch','position',[0 0 12 5])
subplot(121);
for i = -1:1
    for j = -1:1
        viscircles([xCir(:)+i,yCir(:)+j],rCir(:)*gratio);
    end
end
xlim([0 1]); ylim([0 1])
pbaspect([1 1 1]); box on
title('Axon Packing','interpreter','latex','fontsize',20)
set(gca,'xtick',[],'ytick',[])

% Plot the lookup table
subplot(122);
cmap = colormap('parula');
Ibg = A==0;                     % background region
Iol = A>Nax;                    % two-axon region
A2 = ceil(single(A)/Nax*64);    % rescale the colormap for # axons
A2(Ibg) = 1;
A2(Iol) = 1;
[nx,ny] = size(A2);
imgc = cmap(uint16(A2(:)),:);
imgc(Ibg,:) = 0;                % background region is black
imgc(Iol,:) = 1;                % two-axon region is white
imgc = reshape(imgc,[nx,ny,3]);

image(rot90(imgc)); caxis([0 Nax]);
box on; axis off
title('Lookup Table','interpreter','latex','fontsize',20)

%% Run the simulation in the intra-cylindrical space in 2d
% Parameters are defined in the 'simParamInput.txt', and the # particle is
% defined in main.cpp, line 50.
cd(fullfile(root,'cpp'))
system('g++ -std=c++11 RNG.cpp diffusion_lib.cpp main.cpp -o my_code')
system('./my_code')

%% Read simulation results

% read simulation parameters
fileID = fopen(fullfile(root,'cpp','sim_para'),'r');
n=1; C = {};
tline=fgetl(fileID);
while ischar(tline)
    C{n}=tline;
    tline=fgetl(fileID);
    n=n+1;
end
fclose(fileID);

dt = str2double(C{1});        % time step
TN = str2double(C{2});        % # step
Np = str2double(C{3});        % # particles
Tn = str2double(C{4});        % # recorded time step
Din = str2double(C{5});       % IAS diffusivity
Dout = str2double(C{6});      % EAS diffusivity
kappa = str2double(C{7});     % permeability
Dmr = str2double(C{8});       % myelin diffusivity, radial
Dmc = str2double(C{9});       % myelin diffusivity, circumferential
Dmz = str2double(C{10});      % myelin diffusivity, along z-axis
b = str2num(C{end});          % bvalue

% read simulation result
dx2 = load(fullfile(root,'cpp','dx2_diffusion'))/Np;
dy2 = load(fullfile(root,'cpp','dy2_diffusion'))/Np;
dx4 = load(fullfile(root,'cpp','dx4_diffusion'))/Np;
dy4 = load(fullfile(root,'cpp','dy4_diffusion'))/Np;

% calculate radial diffusivity (RD) and radial kurtosis (RK)
t = dt*(1:TN/Tn:TN).';
Dx = dx2/2./t;
Dy = dy2/2./t;
RD = (Dx+Dy)/2;

Kx = dx4./dx2.^2 - 3;
Ky = dy4./dy2.^2 - 3;
RK = (Kx+Ky)/2;

%% Plot time-dependence inside the non-absorbing impermeable cylinder (2d)
%  Particles diffuse in intra-axonal space only.

% inner radius
rin = rCir(1)*gratio*fov;

% Analytical solution of RD and RK inside a cylinder
addpath(fullfile(root,'lib'))
[RD_theory,RK_theory] = cylinderDK(rin,Din,t,200);

figure; set(gcf,'unit','inch','position',[1 1 14 7])
subplot(121)
hold on
h = plot(t,RD/Din,'-');
axis square; grid on
set(h,'linewidth',4)
h_theory = plot(t,RD_theory/Din,'--r');
set(h_theory,'linewidth',3);
legend([h,h_theory],{'Simulation','Exact'},'interpreter','latex','fontsize',30)

xlim([0 max(t)]); ylim([0 1.2]); box on;
set(gca,'fontsize',20,'ytick',0:0.5:2)
xlabel('$t$ (ms)','interpreter','latex','fontsize',30);
ylabel('$D(t)/D_0$','interpreter','latex','fontsize',30);

subplot(122)
hold on
h=plot(t,RK,'-'); axis square
hr = refline(0,0); set(hr,'color','k')
set(h,'linewidth',4,'markerfacecolor',[0 0 1])
h_theory=plot(t,RK_theory,'--r');
set(h_theory,'linewidth',3);

legend([h,h_theory],{'Simulation','Exact'},'interpreter','latex','fontsize',30,'location','southeast')

set(gca,'fontsize',20,'ytick',-2:0.5:2); grid on
xlim([0 sqrt(max(t))]); ylim([-2 0.5]); box on;
xlabel('$t$ (ms)','interpreter','latex','fontsize',30);
ylabel('$K(t)$','interpreter','latex','fontsize',30);
