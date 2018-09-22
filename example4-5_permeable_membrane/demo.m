% This example demonstrates how to calculate the permeability based on the
% particle density and check against the input value, shown in Figure 4,
% point 5 in (Fieremans and Lee, NeuroImage 2018), with more details in
% supplementary information.
%
% Author: Hong-Hsi Lee, September, 2018 (orcid.org/0000-0002-3663-6559)

% ********** Setup the directory on your computer **********
root = 'your_directory_to_this_demo/example4-5_permeable_membrane';
% root = pwd;

% Create files for simualtion
% 1. packing directory
root_packing = fullfile(root,'packing');
fileID = fopen(fullfile(root,'cpp','root.txt'),'w');
fprintf(fileID,root_packing);
fclose(fileID);

% 2. Simulation parameters
dt = 7.5e-3;  % time of each step, ms
TN = 2e4;     % # steps
Din = 2;      % IAS diffusivity, µm^2/ms
Dout = 2;     % EAS diffusivity, µm^2/ms
kappa = 0.05; % permeability, µm/ms
Dmr = 0.2;    % myelin diffusivity, radial, µm^2/ms
Dmc = 0.5;    % myelin diffusivity, circumferential, µm^2/ms
Dmz = 0.5;    % myelin diffusivity, along z-axis, µm^2/ms
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
subplot(121)
hold on
for i = -1:1
    for j = -1:1
        viscircles([xCir(:)+i,yCir(:)+j],rCir(:)*gratio);
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

%% Run the simulation in two coaxial permeable cylinders in 2d
% Parameters are defined in the 'simParamInput.txt', and the # particle is
% defined in main.cpp, line 50. The diffusion is initialized from the
% cylindrical center.
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
xt = load(fullfile(root,'cpp','x_diffusion'));
yt = load(fullfile(root,'cpp','y_diffusion'));
zt = load(fullfile(root,'cpp','z_diffusion'));

% calculate radial diffusivity (RD) and radial kurtosis (RK)
t = dt*(1:TN/Tn:TN).';
Dx = dx2/2./t;
Dy = dy2/2./t;
RD = (Dx+Dy)/2;

Kx = dx4./dx2.^2 - 3;
Ky = dy4./dy2.^2 - 3;
RK = (Kx+Ky)/2;

%% Plot normalized average particle density along radial direction
% Project 2d particle density to 1d particle density

rin = rCir(1)*gratio*fov;   % inner radius
rout = rCir(1)*fov;         % outer radius
xcen = 0.5*fov;             % circle center, x coordinate
ycen = 0.5*fov;             % circle center, y coordinate
nbin = 200;                 % # bin for particle density

% calculate particle density
binRange = linspace(0,fov/2,nbin).';
density = zeros(nbin,Tn);   % particle density
rt = sqrt((xt-xcen).^2+(yt-ycen).^2);
binWidth = binRange(2)-binRange(1);
binCenter = binRange + (binRange(2)-binRange(1))/2;
for i = 1:Tn
    density(:,i) = histc(rt(i,:).',binRange)/Np./(2*pi*binCenter)/binWidth;
end

% plot particle density around outer myelin sheath at t = 90 ms
figure; set(gcf,'unit','inch','position',[1 1 14 7])
subplot(121)
hold on
I90ms = round(90/max(t)*Tn);
h = plot(binCenter-rout,density(:,I90ms),'-');
set(h,'linewidth',3)
xlim([-0.15 0.15]*fov)
ylim([0 0.0075])
hl = plot([0 0],[0 10],'--r'); set(hl,'linewidth',2)
axis square; box on; set(gca,'xtick',(-0.2:0.1:1)*fov,'ytick',0:0.002:0.02,'fontsize',20)
xlabel('$x-x_M$ ($\mu$m)','interpreter','latex','fontsize',30)
ylabel('$n(x)/N_{\rm{par}}$ ($\mu$m$^{-2}$)','interpreter','latex','fontsize',30)
legend(hl,{'membrane'},'interpreter','latex','fontsize',30)

% calculate permeability of outer myelin sheath based on particle density
kappa_out = zeros(Tn,1);    % permeability based on particle density
for i = 1:Tn
    densityi = density(:,i);
    Iout = find(binRange>rout);
    Iout = Iout(1:10);
    bin_out = binRange(Iout);
    
    Imy = find( (binRange>rin).*(binRange<rout) );
    bin_my = binRange(Imy);
    
    density_out = densityi(Iout);
    density_my = densityi(Imy);
    
    A=[bin_out ones(length(density_out),1)];
    tmp=A\density_out;
    gradient_out=tmp(1);

    A=[bin_my ones(length(density_my),1)];
    tmp=A\density_my;
    gradient_my=tmp(1);
    
    myWidth = 10;
    kappa_out(i)=(Dmr*gradient_my)/(mean(density_out)-mean(density_my(end-myWidth+1:end)));
    
end

% plot calculated permeability of outer myelin sheath
subplot(122)
h = plot(t,kappa_out,'-'); set(h,'linewidth',3)
text(t(I90ms),kappa+0.02,'\downarrow','fontsize',50)
xlim([0 150]); ylim([0 0.2]);
set(gca,'xtick',0:50:200,'ytick',0:0.05:1,'fontsize',20)
hr = refline(0,kappa); set(hr,'color','r','linestyle','--','linewidth',2)
axis square
xlabel('$t$ (ms)','interpreter','latex','fontsize',30);
ylabel('$\kappa$ ($\mu$m/ms)','interpreter','latex','fontsize',30);
legend(hr,{'Ideal'},'interpreter','latex','fontsize',30)


