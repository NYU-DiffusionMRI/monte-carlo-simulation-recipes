% Els Fieremans, Hong-Hsi Lee, Physical and numerical phantoms for the
% validation of brain microstructural MRI: A cookbook, NeuroImage 2018
%
% Randomly packed cylinders in Figure 4
%
% The axonal diameter histogram in corpus callosum is based on Fig. 4 in
% (Aboitiz, et al., Brain Research, 1992, 598:143-153)
%
% The fiber packing algorithm is Donev's work in (Donev, et al., J. Comput.
% Phys., 2005, 202:737-764)
%
% ATTENTION: The value of cylindrical volume fraction can be slightly
% different from the input. The actual value needs to be recalcualted based
% on the output packing.

% ********** Setup the directory on your computer **********
root = 'your_directory_to_this_demo/Densely_packed_cylinder';
% root = pwd;
target = fullfile(root,'packing'); mkdir(target)

% input parameters of your packing
maxdensity = 0.76;  % maximal cylindrical volume fraction
                    % cylindrical volume fraction = sum(pi*outer_radius.^2)/fov^2
gratio = 0.585;     % g-ratio, ratio of inner to outer diameter
hardwallBC = 0;     % boundary condition, 0 for periodic, 1 for hard wall
shrinkFactor = 1;   % shrinkage factor, usually > 1
                    % diameter = (diameter in histology)*shrinkfactor

% load inner axonal diameter histogram in corpus callosum
load(fullfile(root,'diameter_histogram','CC_diameter_histogram.mat'));

% use the histogram in genu of the corpus callosum, roi = 1
roi = 1;
N = 500;                                % # axon
pct = cc(roi).frequency;                % frequency, count in percentage
ri = cc(roi).diameter/2*shrinkFactor;   % inner radius
ct = floor(N.*pct);                     % actual count
N = sum(ct);                            % # axon, recalculated
rinit = zeros(N,1);                     % initialize inner radius for packing
i = 0;
for bin = 1:length(ct)
  rinit(i+1:i+ct(bin)) = ri(bin);
  i = i+ct(bin);
end

% plot inner diameter histogram
edges = (0:0.2:9)*shrinkFactor;
Nc = histcounts(2*rinit,edges);
bar(edges(2:end),Nc/N*100,1);
xlim([0 9]); ylim([0 30])
box on; pbaspect([2 1 1])
set(gca,'xtick',0:9,'ytick',0:5:30,'fontsize',12)
title([cc(roi).name,sprintf(', shrinkage factor = %.2f',shrinkFactor)],'fontsize',20)
xlabel('Inner Diameter (µm)','fontsize',16)
ylabel('Frequency (%)','fontsize',16)

%% initialization and generate read.dat for Donev's C++ code
% initial outer radius
rinit = rinit./gratio;

% initial position and rescaled outer radius
[xinit, yinit, rs] = initposition(rinit);

% show initial configuration
figure
viscircles([xinit,yinit],rs,'linewidth',1);
box on; grid on; pbaspect([1 1 1])
xlim([0 1]); ylim([0 1])
set(gca,'xtick',[],'ytick',[])

% prepare input file for C++
C = cell(1,15);
fid = fopen(fullfile(root,'spheres_poly/input'),'r');
for k = 1:15, C{k} = fgetl(fid); end
fclose(fid);

% change N and maxpf variables in inputN file 
C{2}=['int N = ' num2str(N) ';                        // number of spheres'];
C{4}=['double maxpf = ' num2str(maxdensity) ';                  // max packing fraction'];
C{12}=['int hardwallBC = ' num2str(hardwallBC) ';                   // 0 for periodic, 1 for hard wall BC'];
fid = fopen(fullfile(root,'spheres_poly/inputN'), 'w');
fprintf(fid, '%s\n', C{:});
fclose(fid);


%% compile and run Donev's C++ code
cd(fullfile(root,'spheres_poly'))
disp('compiling .............');
system('g++ -O3 -ffast-math -o spheres neighbor.C spheres.C box.C sphere.C event.C heap.C read_input.C');
disp('running C++ ...........');
tic
system('./spheres inputN'); % run
toc
cd ..

%% Save final packing
% read final packing
outputfilename = fullfile(root,'spheres_poly/write.dat');
M = dlmread(outputfilename,'',6,0);
xc = M(:,1); yc = M(:,2); rc = M(:,3)/2;

% calculate fov in µm
res = sum(ri/gratio.*pct)/mean(rc);

% save packing
save(fullfile(target,'packing_parameter.mat'),'xc','yc','rc','gratio','res')

% report cylindrical volume fracion and intra-axonal volume fraction
fprintf('Cylindrical volume fraction = %.4f\n',sum(pi*rc.^2))
fprintf('Intra-axonal volume fraction = %.4f\n',sum(pi*(rc*gratio).^2));

%% Create lookup table
% load packing parameter
load(fullfile(target,'packing_parameter.mat'))

% create lookup table
n = 1000; % size of the lookup table = n x n
[A,B,Nmax] = createlookuptable(n,xc,yc,rc);

% save lookup table
save(fullfile(target,'lookup_table.mat'),'A','n','Nmax')

%% Plot packing and lookup table
load(fullfile(target,'packing_parameter.mat'))
load(fullfile(target,'lookup_table.mat'))

% plot packing
figure; set(gcf,'unit','inch','position',[0 0 12 5])
subplot(121);
for i = -1:1
    for j = -1:1
        viscircles([xc(:)+i,yc(:)+j],rc(:));
    end
end
xlim([0 1]); ylim([0 1])
pbaspect([1 1 1]); box on
title('Axon Packing','interpreter','latex','fontsize',20)
set(gca,'xtick',[],'ytick',[])

% plot lookup table
subplot(122);
cmap = colormap('parula');
Nax = length(rc);               % # axon
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

%% Save files for C++ simulation code
load(fullfile(target,'packing_parameter.mat'))
load(fullfile(target,'lookup_table.mat'))

fid = fopen(fullfile(target,'phantom_res.txt'),'w');
fprintf(fid,'%f',res);
fclose(fid);

fid = fopen(fullfile(target,'phantom_NPix.txt'),'w');
fprintf(fid,'%u',n);
fclose(fid);

fid = fopen(fullfile(target,'phantom_APix.txt'),'w');
for i = 1:size(A,1)
    fprintf(fid,'%u ',A(i,:));
    fprintf(fid,'\n');
end
fclose(fid);

fid = fopen(fullfile(target,'phantom_NAx.txt'),'w');
fprintf(fid,'%u',length(rc));
fclose(fid);

fid = fopen(fullfile(target,'phantom_xCir.txt'),'w');
fprintf(fid,'%f ',xc);
fclose(fid);

fid = fopen(fullfile(target,'phantom_yCir.txt'),'w');
fprintf(fid,'%f ',yc);
fclose(fid);

fid = fopen(fullfile(target,'phantom_rCir.txt'),'w');
fprintf(fid,'%f ',rc);
fclose(fid);

fid = fopen(fullfile(target,'phantom_Nmax.txt'),'w');
fprintf(fid,'%u',Nmax);
fclose(fid);

fid = fopen(fullfile(target,'phantom_gratio.txt'),'w');
fprintf(fid,'%f ',gratio);
fclose(fid);

