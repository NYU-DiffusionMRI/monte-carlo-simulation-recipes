% initialize positions, run Donev's code and show the results

more off; format long; format compact

maxdensity = 0.60;
gratio_min = 0.60;
hardwallBC = 0;   % 0 for periodic, 1 for hard wall boundary conditions
shrink_Fac = 1.5;

if 1
    cd '/Users/magda/Desktop/diffusion_simulation/generate_packing/axon_radius_histogram/aboitiz_human_cc/'
    load cc_hist
    
    pct = cc(1).hi;
    r = cc(1).di/2*shrink_Fac;
    N = 500;
    n = floor(N.*pct);
    N = sum(n);
    rinit = zeros(N,1);
    i = 0;
    Nbins = length(n); 
   for bin=1:Nbins
      rinit(i+1:i+n(bin))=r(bin);
      i=i+n(bin);
   end
    cd ('/Users/magda/Desktop/diffusion_simulation/generate_packing/generate_2d_packings_aboitiz')
end

% if 0
%     cd '/Users/lauren/Documents/simdata/Simbrain/Lamantia distributions/'
%     load hist_sec2
%     
%     pct = perc_sec2;
%     r = d./2;
%     N = 9999;
%     n = N.*pct;
%     rinit = zeros(N,1);
%     i = 0;
%     Nbins = length(n); 
%    for bin=1:Nbins
%       rinit(i+1:i+n(bin))=r(bin);
%       i=i+n(bin);
%    end;
%     cd ('/Users/lauren/Documents/generate_2d_packings')
% end
% 
% if 0
%     cd '/Users/lauren/Documents/simdata/Sciatic'
%     load hist_sciatic
%     
%     pct = perc_sciatic;
%     r = d./2;
%     N = 9999;
%     n = N.*pct;
%     rinit = zeros(N,1);
%     i = 0;
%     Nbins = length(n); 
%    for bin=1:Nbins
%       rinit(i+1:i+n(bin))=r(bin);
%       i=i+n(bin);
%    end;
%     cd ('/Users/lauren/Documents/generate_2d_packings')
% end
    


%% generate read.dat for C++

% % continuous distribution from ISMRM 2012
% if 0 
%    load('ISMRMpacking.mat');
%    N=sum(n); Nbins=length(n);
%    rinit=zeros(N,1);
%    i=0;
%    for bin=1:Nbins
%       rinit(i+1:i+n(bin))=r(bin);
%       i=i+n(bin);
%    end;
% end;
% 
% % a truly bimodal distribution
% if 0
%    N=9999; % total # of disks
%    xi=0.6; % ratio of area of large disks to small disks
%    rS=0.5; rL=2;
%    % xi = NL*rL^2/(NS*rS^2); NL/NS=xi*(rS/rL)^2; NL+NS=N;
%    NL=N*xi*(rS/rL)^2/(1+xi*(rS/rL)^2);
%    NL=round(NL);
%    rinit=rS*ones(N,1);
%    rinit(1:NL)=rL; 
%    Nbins=100; % for further histogram of radii
% end;
% 
% %unimodal distribution
% if 0
%     n = 9999;
%     r = 0.5;
%    N=sum(n); Nbins=length(n);
%    rinit=zeros(N,1);
%    i=0;
%    for bin=1:Nbins
%       rinit(i+1:i+n(bin))=r(bin);
%       i=i+n(bin);
%    end;
% end

rinit = rinit./gratio_min;

[xinit, yinit, rs] = init_positions(rinit);


%% show initial configuration
if 1
   figure; hold on; grid on; axis square; axis equal;
   theta=0:0.01:2*pi; ct=cos(theta); st=sin(theta);
   for i=1:N, plot(xinit(i)+rs(i)*ct,yinit(i)+rs(i)*st,'b-'); end
   
   figure; hist(rinit,Nbins);
end


%% prepare input file for C++
C = cell(1,15);
inpfile0 = fopen('spheres_poly/input','r');
for k=1:15, C{k}=fgetl(inpfile0); end
fclose(inpfile0);
% change N and maxpf variables in inputN file 
C{2}=['int N = ' num2str(N) ';                        // number of spheres'];
C{4}=['double maxpf = ' num2str(maxdensity) ';                  // max packing fraction'];
C{12}=['int hardwallBC = ' num2str(hardwallBC) ';                   // 0 for periodic, 1 for hard wall BC'];
inpfile = fopen('spheres_poly/inputN', 'w');
fprintf(inpfile, '%s\n', C{:});
fclose(inpfile);


%% compile and run C++
cd spheres_poly
disp('compiling .............');
system('g++ -O3 -ffast-math -o spheres neighbor.C spheres.C box.C sphere.C event.C heap.C read_input.C');
disp('running C++ ...........');
tic
system('./spheres inputN'); % run
toc
cd ..


%% show final packing
outputfilename='spheres_poly/write.dat';
Min=dlmread(outputfilename,'',6,0); M=Min;
% M=Min(8:end,:);
x=M(:,1); y=M(:,2); r=M(:,3)/2;
if 1
   figure; hold on; grid on; axis square; axis equal;
   theta=0:0.01:2*pi; ct=cos(theta); st=sin(theta);
   for i=1:N, plot(x(i)+r(i)*ct,y(i)+r(i)*st,'b-'); end
   
   figure; hist(r,Nbins);
end

target = '/Users/magda/Desktop/diffusion_simulation/generate_packing/axon_radius_histogram/aboitiz_human_cc/';
load(fullfile(target,'cc_hist.mat'))
di = cc.di;
ri = di/2*shrink_Fac;
hi = cc.hi;
res = sum(ri/gratio_min.*hi)/sum(hi)/mean(r);

save(['packing_aboitiz_rho' num2str(max_density*100) '_g' num2str(gratio_min*100) '.mat'], 'x', 'y', 'r','gratio_min','res');

%% Pixelize the packing

load(['/Users/magda/Desktop/diffusion_simulation/generate_packing/generate_2d_packings_aboitiz/packing_aboitiz_rho'...
    num2str(max_density*100) '_g' num2str(gratio_min*100) '.mat'])
Np = 1000; % Pixel # along each side
[A,B,NPix,Nmax] = PixelizeGeometry_HHL_v2(Np,x,y,r);
figure; imagesc(A,[0 Nmax]); axis equal; box on; axis off
figure; viscircles([y,-x],r,'linewidth',0.5); axis equal; xlim([0 1]); ylim([-1 0]); box on; axis off

save(['/Users/magda/Desktop/diffusion_simulation/generate_packing/generate_2d_packings_aboitiz/PixGeo_aboitiz_rho'...
    num2str(max_density*100) '_g' num2str(gratio_min*100) '.mat'],'A','NPix','Nmax')

%% Save files for  c++ code

target = ['/Users/magda/Desktop/diffusion_simulation/cpp_practice/packing/packing_aboitiz_rho'...
    num2str(max_density*100) '_g' num2str(gratio_min*100) '_shrinkFac' num2str(shrink_Fac)];
mkdir(target)

fileID = fopen(fullfile(target,'phantom_res.txt'),'w');
fprintf(fileID,'%f',res);
fclose(fileID);

load(['/Users/magda/Desktop/diffusion_simulation/generate_packing/generate_2d_packings_aboitiz/PixGeo_aboitiz_rho'...
    num2str(max_density*100) '_g' num2str(gratio_min*100) '.mat'])
fileID = fopen(fullfile(target,'phantom_NPix.txt'),'w');
fprintf(fileID,'%u',NPix);
fclose(fileID);

fileID = fopen(fullfile(target,'phantom_APix.txt'),'w');
for i = 1:size(A,1)
    fprintf(fileID,'%u ',A(i,:));
    fprintf(fileID,'\n');
end
fclose(fileID);

load(['packing_aboitiz_rho' num2str(max_density*100) '_g' num2str(gratio_min*100) '.mat'])
fileID = fopen(fullfile(target,'phantom_NAx.txt'),'w');
fprintf(fileID,'%u',length(r));
fclose(fileID);

fileID = fopen(fullfile(target,'phantom_xCir.txt'),'w');
fprintf(fileID,'%f ',x);
fclose(fileID);

fileID = fopen(fullfile(target,'phantom_yCir.txt'),'w');
fprintf(fileID,'%f ',y);
fclose(fileID);

fileID = fopen(fullfile(target,'phantom_rCir.txt'),'w');
fprintf(fileID,'%f ',r);
fclose(fileID);

fileID = fopen(fullfile(target,'phantom_Nmax.txt'),'w');
fprintf(fileID,'%u',Nmax);
fclose(fileID);

fileID = fopen(fullfile(target,'phantom_gratio.txt'),'w');
fprintf(fileID,'%f ',gratio_min);
fclose(fileID);

%%

a=[1 2 3;4 5 6];
fileID = fopen(fullfile(target,'test.txt'),'w');
for i = 1:size(a,1)
    fprintf(fileID, '%u ', a(i,:));
    fprintf(fileID, '\n');
end
fclose(fileID);

%%

root='/Users/magda/Desktop/diffusion_simulation/cpp_practice/practice/practice';
fileID = fopen(fullfile(root,'x_diffusion'),'r');
n=1;
tline=fgetl(fileID);
while ischar(tline)
    C{n}=tline;
    tline=fgetl(fileID);
    n=n+1;
end
fclose(fileID);

clear xdiff
for i=1:n-1
    xdiff(i,:)=str2num(C{i});
end

clear tline C
fileID = fopen(fullfile(root,'y_diffusion'),'r');
n=1;
tline=fgetl(fileID);
while ischar(tline)
    C{n}=tline;
    tline=fgetl(fileID);
    n=n+1;
end
fclose(fileID);

clear ydiff
for i=1:n-1
    ydiff(i,:)=str2num(C{i});
end

load('/Users/magda/Desktop/diffusion_simulation/generate_packing/generate_2d_packings_aboitiz/packing_aboitiz_rho76_g585.mat')
close all
figure;
hold on
for i=1:size(xdiff,2)
    plot(xdiff(:,i),ydiff(:,i),'bo-','markersize',5)
    plot(xdiff(1,i),ydiff(1,i),'m.','markersize',10)
end
viscircles([x,y],r,'linewidth',0.5); axis equal; xlim([0 1]); ylim([0 1]); box on;

for i=1:size(xdiff,2)
    m2(i,1)=sum((xdiff(:,i)-xdiff(1,i)).^2+(ydiff(:,i)-ydiff(1,i)).^2);
end

%%

root='/Users/magda/Desktop/diffusion_simulation/cpp_practice/diffusion_myelin_exchange_ellipse_propagator/diffusion_myelin_exchange/';
fileID = fopen(fullfile(root,'x_diffusion'),'r');
n=1;
tline=fgetl(fileID);
while ischar(tline)
    C{n}=tline;
    tline=fgetl(fileID);
    n=n+1;
end
fclose(fileID);

clear xdiff
for i=1:n-1
    xdiff(i,:)=str2num(C{i});
end

clear tline C
fileID = fopen(fullfile(root,'y_diffusion'),'r');
n=1;
tline=fgetl(fileID);
while ischar(tline)
    C{n}=tline;
    tline=fgetl(fileID);
    n=n+1;
end
fclose(fileID);

clear ydiff
for i=1:n-1
    ydiff(i,:)=str2num(C{i});
end

load('/Users/magda/Desktop/diffusion_simulation/generate_packing/generate_2d_packings_aboitiz/packing_aboitiz_rho76_g585.mat')
close all
figure;
hold on
gratio=0.585;
for i=-1:1
    for j=-1:1
        viscircles([x+i,y+j],r,'linewidth',0.5); axis equal; xlim([0 1]); ylim([0 1]); box on;
        viscircles([x+i,y+j],r*gratio,'linewidth',0.5);
    end
end

% axon_label=[371 213];
% viscircles([x(axon_label),y(axon_label)],r(axon_label),'linewidth',0.5); axis equal; xlim([0 1]); ylim([0 1]); box on;
% viscircles([x(axon_label),y(axon_label)],r(axon_label)*gratio,'linewidth',0.5); axis equal; xlim([0 1]); ylim([0 1]); box on;

for i=1:size(xdiff,2)
    plot(xdiff(:,i),ydiff(:,i),'bo-','markersize',5)
    plot(xdiff(1,i),ydiff(1,i),'m.','markersize',10)
    plot(xdiff(end,i),ydiff(end,i),'g.','markersize',10)
end

for i=1:size(xdiff,2)
    m2(i,1)=sum((xdiff(:,i)-xdiff(1,i)).^2+(ydiff(:,i)-ydiff(1,i)).^2);
end


