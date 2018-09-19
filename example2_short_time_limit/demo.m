% Els Fieremans, Hong-Hsi Lee, Physical and numerical phantoms for the
% validation of brain microstructural MRI: A cookbook, NeuroImage 2018
%
% Supplementary Information: Recipes of MC simulations in Figure 4
%
% Example 2: Check short-time limit: Diffusion in extra-cylindrical space
%   of randomly packed impermeable cylinders in 2d

% ********** Setup the directory on your computer **********
root = 'your_directory_to_this_demo/example2_short_time_limit';

% Create files for simualtion
% 1. packing directory
root_packing = fullfile(root,'packing');
fileID = fopen(fullfile(root,'cpp','root.txt'),'w');
fprintf(fileID,root_packing);
fclose(fileID);

% 2. Simulation parameters
dt = 3.5e-5;  % time of each step, ms
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
%   rCir: axon's inner radius
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
        viscircles([xCir(:)+i,yCir(:)+j],rCir(:));
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

%% Run the simulation in the extra-cylindrical space in 2d
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

%% Plot short-time limit for randomly packed impermeable cylinders (2d)
%  Particles diffuse in extra-axonal space only.

figure; set(gcf,'unit','inch','position',[1 1 14 7])
subplot(121)
hold on
h = plot(sqrt(t),RD/Dout,'-'); axis square
set(h,'linewidth',3)

% calculate surface S, volume V, surfave-to-volume ratio SoV
S = fov*2*pi*sum(rCir);
V = fov^2*(1-pi*sum(rCir.^2));
SoV= S/V;
D_theory = Dout*(1-SoV/2*(4/3)*sqrt(Dout*t)/sqrt(pi));
K_theory = (3/8)*(8/5)*SoV*sqrt(Dout*t)/sqrt(pi);

h_theory = plot(sqrt(t),D_theory/Dout,'--r'); set(h_theory,'linewidth',2);
legend([h,h_theory],{'Simulation','Mitra Limit'},'interpreter','latex','fontsize',30,'location','southwest')

set(gca,'xtick',0:0.05:1,'ytick',0:0.5:1,'fontsize',20);
xlim([0 0.3]); 
ylim([0 1.25]); box on; grid on
xlabel('$\sqrt{t}$ (ms$^{1/2}$)','interpreter','latex','fontsize',30);
ylabel('$D(t)/D_0$','interpreter','latex','fontsize',30);

subplot(122);
hold on
h = plot(sqrt(t),RK,'-'); axis square
set(h,'linewidth',3)
h_theory = plot(sqrt(t),K_theory,'--r');
set(h_theory,'linewidth',2);
hl = refline(0,0); set(hl,'color','k')
legend([h,h_theory],{'Simulation','Mitra limit'},'interpreter','latex','fontsize',30,'location','southeast')

set(gca,'xtick',0:0.05:1,'ytick',-2:0.5:1,'fontsize',20);
xlim([0 0.3]); 
ylim([-2 1]); box on; grid on
xlabel('$\sqrt{t}$ (ms$^{1/2}$)','interpreter','latex','fontsize',30);
ylabel('$K(t)$','interpreter','latex','fontsize',30);
